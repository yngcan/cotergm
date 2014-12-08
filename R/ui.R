# FORTESS MONITORING 
# 
# Author: kirk
###############################################################################
library(shiny)
source("chooser.R")

customTextInput<-function (inputId, label, value="",...) {
	tagList(tags$label(label, `for` = inputId), tags$input(id = inputId,
					type="text",
					value=value,...))
}


shinyUI(
		pageWithSidebar(
				
				headerPanel("COSTERGM"),
				sidebarPanel(fluidRow(
						column(3,		h4("net1"),
								radioButtons("net1node1","node1",c("-"=1,"+"=2)),
								radioButtons("net1node2","node2",c("-"=1,"+"=2)),
								radioButtons("net1edge","edge",c("on"=1,"off"=0))),
						column(3,	h4("net2"),
								radioButtons("net2node1","node1",c("-"=1,"+"=2)),
								radioButtons("net2node2","node2",c("-"=1,"+"=2)),
								radioButtons("net2edge","edge",c("on"=1,"off"=0))
					)),
			verbatimTextOutput("test")
				),
				
				mainPanel(
						plotOutput("plot", height = "800px")
				)
		
		
		))
