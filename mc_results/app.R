## app.R ##
library(shiny)
library(shinydashboard)
library(tidyverse)


ui <- dashboardPage(
    
    dashboardHeader(title = "Monte Carlo"),
    dashboardSidebar(
        selectInput(inputId = "prop_missing",
                    label = h5("Select level of missingness (%)"),
                    choices = c(unique(as.numeric(combined_shiny_data$prop_missing)))
        )
    ),
    dashboardBody(box(plotOutput("basic_density")),
                  box(plotOutput("zip_plot")), 
                  box(DT::dataTableOutput('table')), 
                  box(title = "Mechanism of Missingness",width = 2, height = 150,background = "blue",
                      "Missing At Random"), 
                  box(title = "Number of Imputations",width = 2, height = 150, background = "light-blue",
                      "m = 10 for Rubin; m = 2 otherwise"), 
                  box(title = "Number of Subsamples",width = 2, height = 150,background = "aqua",
                      "j = 200")
    )
)

server <- function(input, output) {
    
    # Data Creation
    {results <- 
        read_rds("combined_results_final.rds")
    
    duplicated_names <- duplicated(colnames(results))
    
    results <- results[!duplicated_names]
    
    rubins <- results %>%
        select(.,contains("rub") | prop_missing) %>%
        mutate(type = "Rubin") %>%
        rename(point_estimate = point_estimate_rub, 
               LB = LB_rub, 
               UB = UB_rub, 
               width = rubin_width, 
               prop_missing = prop_missing)
    
    boots <- results %>%
        select(.,contains("boot") | prop_missing) %>%
        mutate(type = "Bootstrap") %>%
        rename(point_estimate = point_estimate_boot, 
               LB = LB_boot, 
               UB = UB_boot, 
               width = boot_width, 
               prop_missing = prop_missing)
    
    jackk <- results %>%
        select(.,-contains(c("boot", "rub")) | prop_missing) %>%
        mutate(type = "Jackknife") %>%
        rename(width = jackknife_width, 
               point_estimate = point_estimate_jackknife, 
               prop_missing = prop_missing)
    
    combined_shiny_data <- bind_rows(rubins, boots, jackk) 
    
    combined_shiny_data$true_var <- rep(2, nrow(combined_shiny_data))}
    
    output$basic_density <- renderPlot({
        reshape2::melt(results[c("point_estimate_rub", "point_estimate_jackknife", "point_estimate_boot", "prop_missing")], id.vars = "prop_missing") %>%
            filter(prop_missing == input$prop_missing) %>%
            ggplot(., 
                   aes(x = variable, fill = variable, color = variable, y = value-2)) + 
            facet_grid(. ~ prop_missing, labeller = labeller(prop_missing = labels)) + 
            ggdist::stat_halfeye(
                alpha = 1,
                adjust = .5, 
                width = .6, 
                .width = 0, 
                justification = -.3, 
                point_colour = NA) + 
            geom_point(
                size = 0.01,
                alpha = 0.05,
                position = position_jitter(
                    seed = 1, width = .1
                )
            ) + 
            geom_boxplot(
                fill = "transparent", 
                color="black",
                width = .25, 
                outlier.shape = NA, 
                show.legend=FALSE
            ) +
            coord_cartesian(ylim = c(-1.7, 1.7), clip = "on") + 
            theme(
                axis.title.x = element_text(size = 16),
                axis.text.x = element_text(size = 14),
                axis.title.y = element_text(size = 16), 
                axis.text.y = element_text(size = 14), 
                axis.ticks.y = element_text(size = 14)) + 
            theme_bw() + 
            theme(axis.title.x = element_blank(), 
                  axis.ticks.x = element_blank(), 
                  axis.text.x = element_blank(),
                  legend.text = element_text(size=14), 
                  legend.title = element_text(size=16),
                  panel.spacing=unit(1.5,"lines"), 
                  strip.text = element_text(
                      size = 14), 
                  plot.title = element_text(size = 20)) + 
            labs(
                y = latex2exp::TeX("$\\widehat{\\beta_1} - \\beta_1"), 
                title = "Density Plots for Point Estimate Bias"
            ) + 
            scale_fill_manual(name = "Method", labels = c("Rubin's Rules", "Jackknife", "Bootstrap"), values = c("#F8766D", "#00BA38", "#619CFF")) +                                                                                                
            scale_color_manual(name = "Method", labels = c("Rubin's Rules", "Jackknife", "Bootstrap"), values = c("#F8766D", "#00BA38", "#619CFF")) +    
            geom_hline(yintercept = 0, linetype="dashed", color = "red") 
        
    })
    
    output$table <- DT::renderDataTable({
        df <- combined_shiny_data %>%
            filter(prop_missing == input$prop_missing) %>%
            mutate("Bias" = point_estimate-true_var) %>%
            mutate("Coverage" = 
                       ifelse(true_var > LB & true_var < UB, TRUE, FALSE)) %>%
            group_by(type) %>%
            summarise("Bias" = round(mean(Bias),3), 
                      "Relative Bias" = round(mean(Bias)/2,3), 
                      "Coverage Probability" = round(sum(Coverage)/{nrow(.)/3},3)) %>%
            rename("Type" = type) %>%
            slice(c(3,2,1))
        DT::datatable(df, options = list(paging = FALSE, searching=FALSE)) %>% 
            DT::formatStyle(1:nrow(df), color = "black")
    })
    
    output$zip_plot <- renderPlot({
        set.seed(234)
        jackk_p <- results %>%
            filter(prop_missing == input$prop_missing) %>%
            sample_n(1e2, replace = FALSE) %>%
            mutate(covers = ifelse(UB > true_var & true_var > LB, "Covers", "Does Not Cover")) %>%
            ggplot(., aes(x = 1:100)) + 
            geom_errorbar(aes(ymin = LB, ymax = UB, color = covers)) + 
            theme_bw() + 
            #coord_flip() + 
            labs(
                x = "", 
                y = "C.I."
            ) + 
            geom_hline(yintercept = results$true_var) + 
            scale_color_manual(name = NULL, values=c("#999999", "#FF0000")) + 
            ggtitle("Jackknife Estimator") + 
            xlim(c(0,100)) + 
            ylim(c(-2.5,5))
        
        set.seed(24123)
        boot_p <- results %>%
            filter(prop_missing == input$prop_missing) %>%
            sample_n(1e2, replace = FALSE) %>%
            mutate(covers = ifelse(UB_boot > true_var & true_var > LB_boot, "Covers", "Does Not Cover")) %>%
            ggplot(., aes(x = 1:100)) + 
            geom_errorbar(aes(ymin = LB_boot, ymax = UB_boot, color = covers)) + 
            theme_bw() + 
            #coord_flip() + 
            labs(
                x = latex2exp::TeX("$i^{th} iteration"), 
                y = "C.I."
            ) + 
            geom_hline(yintercept = results$true_var) + 
            scale_color_manual(name = NULL, values=c("#999999", "#FF0000")) + 
            ggtitle("Bootstrap Estimator") + 
            xlim(c(0,100)) + 
            ylim(c(-2.5,5))
        
        set.seed(234)
        rubin_p <- results %>%
            filter(prop_missing == input$prop_missing) %>%
            sample_n(1e2, replace = FALSE) %>%
            mutate(covers = ifelse(UB_rub > true_var & true_var > LB_rub, "Covers", "Does Not Cover")) %>%
            ggplot(., aes(x = 1:100)) + 
            geom_errorbar(aes(ymin = LB_rub, ymax = UB_rub, color = covers)) + 
            theme_bw() + 
            #coord_flip() + 
            labs(
                x = "", 
                y = "C.I."
            ) + 
            geom_hline(yintercept = results$true_var) + 
            scale_color_manual(name = NULL, values=c("#999999", "#FF0000")) + 
            ggtitle("Rubin's Rules") + 
            xlim(c(0,100)) + 
            ylim(c(-2.5,5))
        
        
        ggpubr::ggarrange(rubin_p, jackk_p, boot_p, nrow = 3, ncol = 1, common.legend = TRUE, legend="bottom")
        
    })
    
}

shinyApp(ui, server)

library(rsconnect)
deployApp()
