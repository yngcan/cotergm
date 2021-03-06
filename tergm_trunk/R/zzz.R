.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("tergm", c("statnet"), FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)
  
  .RegisterMHPs()
  .RegisterConstraintImplications()
}

.RegisterMHPs <- function(){
  ergm.MHP.table("c", "Bernoulli", "atleast",  0, "random", "formationMLE")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd",  0, "random", "formationMLE")
  ergm.MHP.table("c", "Bernoulli", "atleast",  1, "TNT", "formationMLETNT")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd",  1, "TNT", "formationMLETNT")

  ergm.MHP.table("c", "Bernoulli", "atmost",  0, "random", "dissolutionMLE")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd",  0, "random", "dissolutionMLE")
  ergm.MHP.table("c", "Bernoulli", "atmost",  1, "TNT", "dissolutionMLETNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd",  1, "TNT", "dissolutionMLETNT")

  ergm.MHP.table("c", "Bernoulli", "atleast+observed",  0, "random", "formationNonObservedMLE")
  ergm.MHP.table("c", "Bernoulli", "atmost+observed",  0, "random", "dissolutionNonObservedMLE")
  
  #################################

  ergm.MHP.table("c", "Bernoulli", "atleastnonminus",  0, "random", "formationPlusMLE")
  ergm.MHP.table("c", "Bernoulli", "atleastnonplus",  0, "random", "formationMinusMLE")
  ergm.MHP.table("c", "Bernoulli", "atmostnonminus",  0, "random", "dissolutionPlusMLE")
  ergm.MHP.table("c", "Bernoulli", "atmostnonplus",  0, "random", "dissolutionMinusMLE")
  
  #################################
  
  ergm.MHP.table("c", "Bernoulli", "atleast+blockdiag",  0, "random", "formationMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd+blockdiag",  0, "random", "formationMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atleast+blockdiag",  1, "TNT", "formationMLEblockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd+blockdiag",  1, "TNT", "formationMLEblockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+blockdiag",  0, "random", "dissolutionMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd+blockdiag",  0, "random", "dissolutionMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atmost+blockdiag",  1, "TNT", "dissolutionMLEblockdiagTNT")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd+blockdiag",  1, "TNT", "dissolutionMLEblockdiagTNT")

  ergm.MHP.table("c", "Bernoulli", "atleast+blockdiag+observed",  0, "random", "formationNonObservedMLEblockdiag")
  ergm.MHP.table("c", "Bernoulli", "atmost+blockdiag+observed",  0, "random", "dissolutionNonObservedMLEblockdiag")

  ergm.MHP.table("f", "Bernoulli", "",  0, "random", "formation")
  ergm.MHP.table("f", "Bernoulli", "bd",  0, "random", "formation")
  ergm.MHP.table("f", "Bernoulli", "",  1, "TNT", "formationTNT")
  ergm.MHP.table("f", "Bernoulli", "bd",  1, "TNT", "formationTNT")
  ergm.MHP.table("d", "Bernoulli", "",  0, "random", "dissolution")
  ergm.MHP.table("d", "Bernoulli", "bd",  0, "random", "dissolution")
  ergm.MHP.table("d", "Bernoulli", "",  1, "TNT", "dissolutionTNT")
  ergm.MHP.table("d", "Bernoulli", "bd",  1, "TNT", "dissolutionTNT")

}

.RegisterConstraintImplications <- function(){
  ergm.ConstraintImplications("atleast", c())
  ergm.ConstraintImplications("atmost", c())
  ergm.ConstraintImplications("nonminus", c())
  ergm.ConstraintImplications("nonplus", c())
}
