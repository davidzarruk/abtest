
#----------------------------------------#
#------------    Funciones    -----------#
#----------------------------------------#

library("lmtest")
library("sandwich")
library('tidyverse')

options(scipen=10000000)

resultados = function(results, CI, obs, req_observaciones, treatments, treatment_variable, control_group_label, factor_mult = 100){
  
  num_treats = length(treatments)
  
  res = data.frame(treatment = array(0, dim=c(num_treats+1)), 
                   teffect = array(0, dim=c(num_treats+1)), 
                   pvalues = array(0, dim=c(num_treats+1)), 
                   sd = array(0, dim=c(num_treats+1)), 
                   CIl = array(0, dim=c(num_treats+1)), 
                   CIh = array(0, dim=c(num_treats+1)), 
                   obs = array(0, dim=c(num_treats+1)), 
                   req_obs = array(0, dim=c(num_treats+1)))
  
  for(i in 1:num_treats){
    treat = treatments[i]
    rowname = paste0("factor(", treatment_variable, ")", treat)
    
    res$teffect[i+1]   = round(factor_mult*results[rowname,1], 3)
    res$pvalues[i+1]   = round(results[rowname,4], 3)
    res$treatment[i+1] = treat
    res$sd[i+1]        = round(results[rowname,2], 3)
    res$CIl[i+1]       = round(factor_mult*CI[rowname,1], 3)
    res$CIh[i+1]       = round(factor_mult*CI[rowname,2], 3)
    res$obs[i+1]       = obs[1, as.character(treat)]
    res$req_obs[i+1]   = req_observaciones[1, as.character(treat)]
  }

  res$teffect[1]   = round(factor_mult*results["(Intercept)",1], 3)
  res$pvalues[1]   = "-"
  res$treatment[1] = control_group_label
  res$sd[1]        = "-"
  res$CIl[1]       = "-"
  res$CIh[1]       = "-"
  res$obs[1]       = obs[1, control_group_label]
  res$req_obs[1]   = "-"
  
  res[2:(num_treats + 1),] = res[2:(num_treats + 1),][order(res[2:(num_treats + 1),]$treatment),]
  
  return(res)
}

grafico_res = function(num_treats, result_data, obs, labels, titulo, axistitle, n_power, power_constraint, add_legend = FALSE){
  
  result_data$color = ""
  if(power_constraint == TRUE){
    try(result_data[as.double(result_data$pvalues) <= 0.05 & as.double(gsub(",", "", result_data$obs)) < as.double(gsub(",", "", result_data$req_obs)),]$color <- "significative", silent=TRUE)
    try(result_data[as.double(result_data$pvalues) <= 0.05 & as.double(gsub(",", "", result_data$obs)) >= as.double(gsub(",", "", result_data$req_obs)),]$color <- "significative_conf", silent=TRUE)
    try(result_data[as.double(result_data$CIl) <= 0 & as.double(result_data$CIh) >= 0 & as.double(gsub(",", "", result_data$obs)) >= n_power,]$color <- "non_sign", silent=TRUE)
    try(result_data[as.double(result_data$CIl) <= 0 & as.double(result_data$CIh) >= 0  & as.double(gsub(",", "", result_data$obs)) < n_power,]$color <- "non_sign_power", silent=TRUE)
  } else{
    try(result_data[as.double(result_data$pvalues) <= 0.05,]$color <- "significative_conf", silent=TRUE)
    try(result_data[as.double(result_data$CIl) <= 0 & as.double(result_data$CIh) >= 0,]$color <- "non_sign", silent=TRUE)
  }
  
  graph = ggplot(result_data, 
                 aes(x = factor(treatment), 
                     y = teffect, 
                     fill = color)) + 
    geom_bar(stat="identity", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=(as.double(CIl)), ymax=(as.double(CIh))), width=.2,
                  position=position_dodge(.9)) + 
    theme_minimal() + 
    geom_hline(yintercept = 0, color = 'black', size = 1) + 
    xlab("") + ylab(axistitle) + theme(plot.title = element_text(hjust=0.5)) + 
    scale_x_discrete(breaks=result_data$treatment,
                     labels=labels) + ggtitle(titulo) +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))
  
  if(add_legend == FALSE){
    graph = graph + theme(legend.position = 'none') 
  }

  if(power_constraint == FALSE){
    graph = graph +
      scale_fill_manual(name = "Significance:",
                        values=c("significative_conf" = "red", 
                                 "non_sign" = "grey50"),
                        labels=c("significative_conf" = "Significant", 
                                 "non_sign" = "Not significant"))  
  } else{
    graph = graph +  
      scale_fill_manual(name = "Significance:",
                        values=c("significative" = "orange", 
                                 "significative_conf" = "red", 
                                 "non_sign_power" = "lightblue", 
                                 "non_sign" = "grey50"),
                        labels=c("significative" = "Significant - low precision", 
                                 "significative_conf" = "Significant - high precision", 
                                 "non_sign_power" = "Not significant - more obs required", 
                                 "non_sign" = "Not significant"))  
  }
  
  return(graph)
}

tabla_res = function(num_treats, result_data, df){
  names(result_data) = c("Tratamiento", "Efecto", "p-valor", "Desv. Est", "IC-botom", "IC-top", "Obs.")
  
  g = tableGrob(result_data, rows = NULL)
  g = gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 2, b = nrow(g), l = 1, r = ncol(g))
  g = gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, l = 1, r = ncol(g)) 
  
  return(g)  

}

confint.robust <- function (object, parm, level = 0.95, ...){
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- stats:::format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,
                                                             pct))
  ses <- sqrt(diag(sandwich::vcovHC(object, type="HC1")))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}

regression_table = function(dataframe, 
                            titulo, 
                            dependent, 
                            treatment_variable = "treatment",
                            controls = c(), 
                            treatment_labels, 
                            control_group_label,
                            factor_mult = 100, 
                            subtitulos = c(), 
                            axistitle = "Puntos porcentuales",
                            significance_level = 0.95,
                            power_constraint = FALSE,
                            power_value = NA,
                            power_delta = NA,
                            power_sd = NA,
                            tolerance_level = 0.0005,
                            tolerance_perc = 0.95,
                            robust_se = TRUE){
  
  if(power_constraint == TRUE){
    n_power = power.t.test(power=power_value, 
                           delta=power_delta,
                           sd=power_sd,
                           type="two.sample",
                           sig.level = significance_level)$n
  } else{
    n_power = 0
  }
  
  if( any(treatment_labels == '') ) stop('The treatment has empty values.')
  
  control_group_label = as.character(control_group_label)
  
  treatments = setdiff(treatment_labels, control_group_label)
  num_treats = length(treatments)
  
  # Output
  output = data.frame(treatment = array(0, dim=c(num_treats+1)), teffect = array(0, dim=c(num_treats+1)), 
                      pvalues = array(0, dim=c(num_treats+1)), sd = array(0, dim=c(num_treats+1)), 
                      CIl = array(0, dim=c(num_treats+1)), CIh = array(0, dim=c(num_treats+1)), obs = array(0, dim=c(num_treats+1)))
  
  labels = array("", dim = c(num_treats))
  
  data = dataframe %>% select(c(dependent, treatment_variable, controls))
  data[[treatment_variable]] = factor(data[[treatment_variable]])
  data = na.omit(data)
  
  # Set base to control variable
  data[[treatment_variable]] = relevel(data[[treatment_variable]], ref = control_group_label)
  
  # Regresion: outcome vs treatment y controles
  
  controls_factors = c("day", "prime", "dapa", "bapa", "platform", "segment", "churn", "maturity", "gama")
  
  controls_with_factors = controls
  for(j in 1:length(controls_factors)){
    if(sum(controls_factors[j] == controls_with_factors) == 1){
      controls_with_factors[controls_factors[j] == controls_with_factors] = paste0("factor(", controls_factors[j],")")
    }
  }
  
  if(length(controls_with_factors)>0){
    model = lm(as.formula(paste(dependent, "~ factor(", treatment_variable, ") + ", paste(controls_with_factors, collapse=" + "))), data = data) 
  } else{
    model = lm(as.formula(paste(dependent, "~ factor(", treatment_variable, ")")), data = data) 
  }

  results = coeftest(model, vcov = vcovHC(model, type="HC1"))
  
  if(robust_se == TRUE){
    CI = confint.robust(model, level = significance_level)
  } else{
    CI = confint(model, level = significance_level)
  }

  observaciones = array(0, dim=c(1,(num_treats+1)))  
  colnames(observaciones) = c(as.character(control_group_label), as.character(sort(treatments)))

  req_observaciones = array(0, dim=c(1,(num_treats)))  
  colnames(req_observaciones) = as.character(sort(treatments))
  
  labels = array("", dim=c(1,num_treats))  
  colnames(labels) = sort(treatments)
  
  for(j in 1:num_treats){
    observaciones[1,as.character(treatments[j])] = format(sum(data[[treatment_variable]] == treatments[j]), big.mark=",")

    # Labels for graph
    if(length(subtitulos) > 0){
      labels[1, as.character(treatments[j])] = paste(subtitulos[j], "\n\nN=", observaciones[1,as.character(treatments[j])])
    } else{
      labels[1, as.character(treatments[j])] = paste(treatments[j], "\n\nN=", observaciones[1,as.character(treatments[j])])
    }
    
    # Observaciones requeridas para estimador confiable
    s = sd(data[data[[treatment_variable]] == treatments[j], ][[dependent]])
    req_observaciones[1, as.character(treatments[j])] = format(round(((1.96*s)/tolerance_level)^2), big.mark=",")
    
  }
  observaciones[1,control_group_label] = format(sum(data[[treatment_variable]] == control_group_label), big.mark=",")

  output = resultados(results, CI, observaciones, req_observaciones, treatments, treatment_variable, control_group_label, factor_mult)

  ggraph = grafico_res(num_treats, output[2:(num_treats+1),], observaciones, labels, titulo, axistitle, n_power, power_constraint)
  
  return(list(output, ggraph))
}
