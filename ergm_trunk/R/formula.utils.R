###################################################################
## This file has utilities whose primary purpose is examining or ##
## manipulating ERGM formulas.                                   ##
###################################################################

## This function appends a list of terms to the RHS of a formula. If
## the formula is one-sided, the RHS becomes the LHS, if
## keep.onesided==FALSE (the default).
## For example,
## append.rhs.formula(y~x,list(as.name("z1"),as.name("z2"))) -> y~x+z1+z2
## append.rhs.formula(~y,list(as.name("z"))) -> y~z
## append.rhs.formula(~y+x,list(as.name("z"))) -> y+x~z
## append.rhs.formula(~y,list(as.name("z")),TRUE) -> ~y+z
append.rhs.formula<-function(object,newterms,keep.onesided=FALSE){
  for(newterm in newterms){
    if(length(object)==3) object[[3]]<-call("+",object[[3]],newterm)
    else if(keep.onesided) object[[2]]<-call("+",object[[2]],newterm)
    else object[[3]]<-newterm
  }
  object
}

# A reimplementation of update.formula() that does not simplify.  Note
# that the resulting formula's environment is set as follows. If
# from.new==FALSE, it is set to that of object. Otherwise, a new
# sub-environment of object, containing, in addition, variables in new
# listed in from.new (if a character vector) or all of new (if TRUE).
ergm.update.formula<-function (object, new, ..., from.new=FALSE){
  old.lhs <- if(length(object)==2) NULL else object[[2]]
  old.rhs <- if(length(object)==2) object[[2]] else object[[3]]
  
  new.lhs <- if(length(new)==2) NULL else new[[2]]
  new.rhs <- if(length(new)==2) new[[2]] else new[[3]]
  
  sub.dot <- function(c, dot){
    if(is.null(dot)) c # If nothing to substitute with, just return it.
    else if(is.call(c)) as.call(c(list(c[[1]]), lapply(c[-1], sub.dot, dot))) # If it's a call, construct a call consisting of the call and each of the arguments with the substitution performed, recursively.
    else if(is.name(c) && c==".")  dot # If it's a dot, return substitute.
    else c # If it's anything else, just return it.
  }
  
  deparen<- function(c, ops = c("+","*")){
    if(is.call(c)){
      if(as.character(c[[1]]) %in% ops){
        op <- as.character(c[[1]])
        if(length(c)==2 && is.call(c[[2]]) && c[[2]][[1]]==op)
          return(deparen(c[[2]], ops))
        else if(length(c)==3 && is.call(c[[3]]) && c[[3]][[1]]==op)
          return(call(op, call(op, deparen(c[[2]],ops), deparen(c[[3]][[2]],ops)), deparen(c[[3]][[3]],ops)))
      }
      return(as.call(c(list(c[[1]]), lapply(c[-1], deparen, ops)))) # If it's a non-reducible call, construct a call consisting of the call and each of the arguments with the substitution performed, recursively.
    }else return(c)
  }
  
  out <- if(length(new)==2) call("~", deparen(sub.dot(new.rhs, old.rhs))) else call("~", deparen(sub.dot(new.lhs, old.lhs)), deparen(sub.dot(new.rhs, old.rhs)))

  #  a new sub-environment for the formula, containing both
  # the variables from the old formula and the new.
  
  if(identical(from.new,FALSE)){ # The new formula will use the environment of the original formula (the default).
    e <- environment(object)
  }else{
    # Create a sub-environment also containing variables from environment of new.
    e <- new.env(parent=environment(object))
    
    if(identical(from.new,TRUE)) from.new <- ls(pos=environment(new)) # If TRUE, copy all of them (dangerous!).
    
    for(name in from.new)
      assign(name, get(name, pos=environment(new)), pos=e)
  }

  as.formula(out, env = e)
}

term.list.formula<-function(rhs, sign=+1){
  if(length(rhs)==1) {attr(rhs,"sign")<-sign; list(rhs)}
  else if(length(rhs)==2 && rhs[[1]]=="+") term.list.formula(rhs[[2]],sign)
  else if(length(rhs)==2 && rhs[[1]]=="-") term.list.formula(rhs[[2]],-sign)
  else if(length(rhs)==3 && rhs[[1]]=="+") c(term.list.formula(rhs[[2]],sign),term.list.formula(rhs[[3]],sign))
  else if(length(rhs)==3 && rhs[[1]]=="-") c(term.list.formula(rhs[[2]],sign),term.list.formula(rhs[[3]],-sign))
  else if(rhs[[1]]=="(") term.list.formula(rhs[[2]], sign)
  else {attr(rhs,"sign")<-sign; list(rhs)}
}


copy.named<-function(x){
  y<-list()
  for(name in names(x)) y[[name]]<-x[[name]]
  y
}


model.transform.formula <- function(object, theta, response=NULL, recipes, ...){
  ## Recipe syntax:
  ##
  ## Recipes are a named list with the names representing a map from
  ## term name to recipe. Each recipe is a list with instructions
  ## about what should be done with a term with that name.
  ##
  #### General settings
  ##
  ## filter: A function that takes a list containing the term name and
  ## arguments to the term, and returns either TRUE or FALSE. If it
  ## returns FALSE, the term is not processed further and is treated
  ## like a term not in the recipes list. If the function is absent,
  ## TRUE is assumed (i.e., all terms with this name are processed).
  ##
  ## custom: A function that takes two arguments: a list comprising
  ## the term name and arguments to the term and a vector of model
  ## parameters for the term. It must return a list with named
  ## elements: `theta`, a vector of model parameters for the term to
  ## be added to the ouptut parameters and `term`, an unevaluated call
  ## (or a list with the term name being the first element and
  ## arguments being the remaining (named) elements) to be appended to
  ## the rest of the formula. *If present, overrides all settings
  ## listed below. Most of the time, you should use those.*
  ##
  #### Convenience settings
  ##
  ## name: A string or a function. If a string, the term will be
  ## renamed to this. If a function, it must take two arguments: a
  ## string with the current name of the term and a list with the term name
  ## and arguments to the term, and return a string with the new name.
  ##
  ## tocoef: Either a numeric vector of integers or a function. If
  ## numeric, gives the indices of the model parameters for the term
  ## that are to be copied directly the output parameter vector. If a
  ## function, the function must take a vector of model parameters for
  ## the term and a list with the term name and arguments to the term
  ## and return a numeric vector to be appended to the output
  ## parameters.
  ##
  ## toarg: A named list of either numeric vectors of integers or
  ## functions (may be heterogeneous). For each element in the list,
  ## the term's corresponding argument will be set to either the model
  ## parameters indicated by the indices (if numeric) or the return
  ## value (if function). As with the others, the function must take
  ## two arguments: a vector of model parameters for the term and a
  ## list with the term name and arguments to the term.
  ##
  ## constant: A named list of elements of any type. The corresponding
  ## arguments of the term are set to the values in the list. This is
  ## a simple special case of toarg, if it were given a function that
  ## returned a constant value.

  m <- ergm.getmodel(object, ergm.getnetwork(object), response=response)
  theta.inds<-cumsum(c(1,coef.sublength.model(m)))
  terms<-term.list.formula(object[[3]])
  form<-object
  ## This deletes the formula's RHS, and LHS becomes RHS (for the moment).
  form[[3]]<-NULL
  newtheta<-c()
  for(i in seq_along(terms)){
    if(!is.call(terms[[i]]) ||
       !(as.character(terms[[i]][[1]]) %in% names(recipes)) ||
       (!is.null(recipes[[as.character(terms[[i]][[1]])]]$filter) &&
        !recipes[[as.character(terms[[i]][[1]])]]$filter(as.list(terms[[i]])))){
      ## If it's not a call OR is a call but does not have a recipe OR
      ## does have a recipe, but the filter function says it should be
      ## skipped (e.g. it's already fixed), then just append it.
      form<-append.rhs.formula(form,list(terms[[i]]))
      newtheta<-c(newtheta,theta[theta.inds[i]:(theta.inds[i+1]-1)])
    }else{
      ## Otherwise, it gets complicated...
      recipe<-recipes[[as.character(terms[[i]][[1]])]]
      orig.list<-call.list<-as.list(terms[[i]])

      if(!is.null(recipe$custom)){
        ## Custom recipe
        out<-recipe$custom(orig.list,theta[theta.inds[i]:(theta.inds[i+1]-1)])
        newtheta<-c(newtheta,out$theta)
        form<-append.rhs.formula(form,
                                 if(is.call(out$term)) list(out$term)
                                 else as.call(c(as.name(out$term[[1]]),out$term[-1])))
      }else{

        if("" %in% names(call.list)[-1]) stop("Curved terms must have all their arguments passed by name.")
        
        ## Rename the term.
        if(!is.null(recipe$name)){
          call.list[[1]]<-
            if(is.function(recipe$name)) as.name(recipe$name(as.character(orig.list[[1]]),orig.list))
            else as.name(recipe$name)
        }
        ## Now, go through the arguments to be replaced:
        ## The constants:
        for(name in names(recipe$constant))
          call.list[[name]]<-recipe$constant[[name]]
        
        ## The elements of theta:
        for(name in names(recipe$toarg))
          call.list[[name]]<-
            if(is.function(recipe$toarg[[name]])) recipe$toarg[[name]](theta[theta.inds[i]:(theta.inds[i+1]-1)],orig.list)
            else theta[theta.inds[i]+recipe$toarg[[name]]-1]
        
        ## Now, add the newly rewritten call to the formula.
        form<-append.rhs.formula(form,list(as.call(call.list)))
        
        ## The parts that remain in theta:
        newtheta<-c(newtheta,
                    if(is.function(recipe$tocoef)) recipe$tocoef[[name]](theta[theta.inds[i]:(theta.inds[i+1]-1)],orig.list)
                    else theta[theta.inds[i]+recipe$tocoef-1])
      }
    }
  }
  list(formula=form,theta=newtheta)
}


## Convert a fitted curved ERGM (or its formula + theta) into a linear
## (fixed=TRUE) model with the curved parameters in theta substituted
## into the formula according to a set of recipes. Returns the new
## formula and the appropriate parameter vector.

fix.curved <- function(object, ...) UseMethod("fix.curved")

fix.curved.ergm <- function(object,...){
  fix.curved.formula(object$formula, coef(object), response=object$response, ...)
}

fix.curved.formula <- function(object, theta, response=NULL, ...){
  recipes<-list()
  is.fixed.1<-function(a) is.null(a$fixed) || a$fixed==FALSE
  recipes$gwdsp<-recipes$gwesp<-recipes$gwnsp<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(alpha=2), constant=list(fixed=TRUE))
  recipes$altkstar<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(lambda=2), constant=list(fixed=TRUE))
  recipes$gwb1degree<-recipes$gwb2degree<-recipes$gwdegree<-recipes$gwidegree<-recipes$gwodegree<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(decay=2), constant=list(fixed=TRUE))

  model.transform.formula(object, theta, response=response, recipes, ...)
}


## Convert a fitted curved ERGM (or its formula + theta) into a curved
## model suitable for use as input to ergm().  This is a workaround
## around a current issue in ERGM and may be eliminated in the future.

enformulate.curved <- function(object, ...) UseMethod("enformulate.curved")

enformulate.curved.ergm <- function(object,...){
  fix.curved.formula(object$formula, coef(object), response=object$response, ...)
}

enformulate.curved.formula <- function(object, theta, response=NULL, ...){
  recipes<-list()
  is.fixed.1<-function(a) is.null(a$fixed) || a$fixed==FALSE
  recipes$gwdsp<-recipes$gwesp<-recipes$gwnsp<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(alpha=2))
  recipes$altkstar<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(lambda=2))
  recipes$gwb1degree<-recipes$gwb2degree<-recipes$gwdegree<-recipes$gwidegree<-recipes$gwodegree<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(decay=2))

  model.transform.formula(object, theta, response=response, recipes, ...) 
}

set.offset.formula <- function(object, which, response=NULL){
  nw <- ergm.getnetwork(object)
  m<-ergm.getmodel(object, nw, response=response,role="target")
  to_offset <-unique(rep(seq_along(m$terms),coef.sublength.model(m))[which]) # Figure out which terms correspond to the coefficients to be offset.
  terms <- term.list.formula(object[[3]])
  for(i in to_offset)
    if(!inherits(terms[[i]],"call") || terms[[i]][[1]]!="offset") # Don't offset terms already offset.
      terms[[i]]<-call("offset", terms[[i]]) # Enclose the term in an offset.
  ergm.update.formula(object, append.rhs.formula(~.,terms)) # append.rhs.formula call returns a formula of the form .~terms[[1]] + terms[[2]], etc.
}

unset.offset.formula <- function(object, which=TRUE, response=NULL){
  nw <- ergm.getnetwork(object)
  m<-ergm.getmodel(object, nw, response=response,role="target")
  to_unoffset <-unique(rep(seq_along(m$terms),coef.sublength.model(m))[which]) # Figure out which terms correspond to the coefficients to be un offset.
  terms <- term.list.formula(object[[3]])
  for(i in to_unoffset)
    if(inherits(terms[[i]],"call") && terms[[i]][[1]]=="offset") # Is the term an offset?
      terms[[i]]<-terms[[i]][[2]] # Grab the term inside the offset.
  ergm.update.formula(object, append.rhs.formula(~.,terms)) # append.rhs.formula call returns a formula of the form .~terms[[1]] + terms[[2]], etc.
}

# Delete all offset() terms in an ERGM formula.
remove.offset.formula <- function(object, response=NULL){
  terms <- term.list.formula(object[[3]])
  for(i in rev(seq_along(terms)))
    if(inherits(terms[[i]],"call") && terms[[i]][[1]]=="offset") # Is the term an offset?
      terms[[i]]<-NULL # Delete the offset term.
  ergm.update.formula(object, append.rhs.formula(~.,terms)) # append.rhs.formula call returns a formula of the form .~terms[[1]] + terms[[2]], etc.
}

# A lightweight function that simply returns the offset vectors
# associated with a formula.
offset.info.formula <- function(object, response=NULL){
  nw <- ergm.getnetwork(object)
  m<-ergm.getmodel(object, nw, response=response,role="target")
  with(m$etamap, list(term=offset, theta=offsettheta,eta=offsetmap))
}
