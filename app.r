library(corpcor);library(Matrix);library(shiny)
eigen_dim<-function(mat){
  e.vals<-eigen(mat)$values
  e.vals[which(e.vals<0)]<-0
  sum(e.vals/e.vals[1])
}
eigen_out<-function(mat){
  e.vals<-eigen(mat)$values
  e.vals[which(e.vals<0)]<-0
  c(e.vals[1],e.vals[2],e.vals[2]/e.vals[1],sum(e.vals/e.vals[1]))
}

ui <- fluidPage(
  titlePanel("Calculate Dimensionality Metrics"),
  
  textAreaInput(inputId = "Mat",label="Paste matrix here",width="100%",height="50%"),
  numericInput(inputId = "k",label="Number of traits",value=3,width = "25%"),
  radioButtons(inputId="by_row",label="By Row?",
               c("Yes"="T",
                 "No"="F")),
  actionButton("go", "Import!"),
  tableOutput(outputId = "Orig_Mat"),
  radioButtons(inputId="triangle",label="Upper or Lower Triangle",
               c("Upper"="U",
                 "Lower"="L")),
  actionButton("do", "Convert!"),
  tableOutput(outputId="New_Mat"),
  actionButton("now", "Calculate!"),
  tableOutput(outputId = "Metrics")
)

server <- function(input, output){
  values <- reactiveValues(mat2 = NULL)
  values <- reactiveValues(mat3=NULL)
  
  observeEvent(input$go,{
    output$Orig_Mat <- renderTable({
      values$mat2 <- matrix(scan(text=input$Mat),input$k,input$k,byrow=input$by_row)
      values$mat2},digits=3)})
  
  observeEvent(input$do, {
    output$New_Mat <- renderTable({
      values$mat3 <- forceSymmetric(values$mat2,input$triangle)
      as.matrix(values$mat3)},digits=3)    })  
  
  observeEvent(input$now,{
    output$Metrics <- renderTable({
      rbind(c(eigen_out(as.matrix(values$mat3))[1:2],
              eigen_dim(as.matrix(values$mat3))))
  },digits=3)})
  
  
}

shinyApp(ui=ui, server=server)
