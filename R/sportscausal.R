sportscausal <-
  function(y.exp, y.con, pre.period, post.period, 
           is.plot = TRUE, model.select = "AIC",
           max.p = 3, max.d = 3, max.q = 3, feature = NULL){
    
    oldpar <- graphics::par(no.readonly = TRUE)   ## to make sure that the par environment does not change after function is called
    on.exit(graphics::par(oldpar)) 
    
    if(length(y.con) == 0){
      stop("empty vector of response in control group")
    }
    
    if(!(length(y.con) == (length(pre.period) + length(post.period)))){
      stop("length of input vector is different")
    }
    
    if(!model.select == "lstm"){
      
      ## input quality control
      
      if(((max.p +1)*(max.d + 1)*(max.q + 1) >= 100) & (model.select == "AIC")){
        warning("too many ARIMA candidates for model selection")
      }
      
      if((model.select == "CV") & (length(y.con) < 50)){
        warning("Cross validation is not suggested to select time series model with small number of data points")
      }
      
      
      if(model.select == "AIC"){
        
        mm = NULL
        
        ## I think the for-loop would be fine since the calculation is light
        
        for (p in 0:max.p) {
          for (d in 0:max.d) {
            for (q in 0:max.q) {
              
              fit.pdq = suppressWarnings(stats::arima(y.con[pre.period], order = c(p,d,q)))
              
              aic.pdq = stats::AIC(fit.pdq)
              
              mm = rbind(mm, c(p,d,q,aic.pdq))
            }
          }
        }
        
        colnames(mm) = c("p", "d", "q", "AIC")
        
        oo = which.min(mm[,4]) ## select the smallest AIC
        
        fit.best = suppressWarnings(stats::arima(y.con[pre.period], order = c(mm[oo,1],
                                                                              mm[oo,2],
                                                                              mm[oo,3])))
        
        predict.best = stats::predict(fit.best, length(post.period))
        
        y0 = c(y.con[pre.period], as.numeric(predict.best$pred))
      }
      
      if(model.select == "CV"){
        
        ## I will use y.con[pre.period] data to train the model, using CV
        
        xx = y.con[pre.period]
        
        ## I split the data into training data set and validation set, 70-30
        
        ss = round(0.7*length(xx))
        
        xx.train = xx[1:ss]
        xx.test = xx[(ss+1):length(xx)]
        
        mm = NULL
        
        ## I think the for-loop would be fine since the calculation is light
        
        for (p in 0:max.p) {
          for (d in 0:max.d) {
            for (q in 0:max.q) {
              
              fit.pdq = suppressWarnings(stats::arima(xx.train, order = c(p,d,q)))
              predict.pdq = stats::predict(fit.pdq, length(xx.test))
              mse.pdq = sum((xx.test - predict.pdq$pred)^2)
              mm = rbind(mm, c(p,d,q,mse.pdq))
            }
          }
        }
        
        colnames(mm) = c("p", "d", "q", "MSE")
        
        oo = which.min(mm[,4]) ## select the smallest MSE
        
        fit.best = suppressWarnings(stats::arima(y.con[pre.period], order = c(mm[oo,1],
                                                                              mm[oo,2],
                                                                              mm[oo,3])))
        
        predict.best = stats::predict(fit.best, length(post.period))
        
        y0 = c(y.con[pre.period], as.numeric(predict.best$pred))
        
      }
      
    }
    
    if(model.select == "lstm"){
      
      if(!is.null(feature)){
        
        ## input quality control
        
        if(!(length(y.con) == dim(feature)[1])){
          stop("length of input vector is different")
        }
        
        if(length(y.con) <= 200){
          stop("sample size is small, lstm prediction might be problematic")
        }
        
        ## make data a copy
        y.con0 = y.con
        feature0 = feature
        
        ## for time series forcasting, I will only need the data from pre period
        y.con = y.con0[pre.period]
        feature = as.matrix(feature0[pre.period,])
        
        ## see reference http://datasideoflife.com/?p=1171
        
        #################################################################
        ####################### data processing #########################
        #################################################################
        
        ## step 1: calculate the scale factors, for the response y.con and also feature
        ## also, normalize the response and feature
        
        scale.factor.response = c(mean(y.con), stats::sd(y.con))
        y.con.norm = (y.con - scale.factor.response[1])/scale.factor.response[2]
        
        scale.factor.feature = matrix(NA, 2, dim(feature)[2])
        feature.norm = feature
        
        for (i in 1:dim(feature)[2]) {
          
          xx.i = feature[,i]
          mu.i = mean(xx.i)
          sd.i = stats::sd(xx.i)
          scale.factor.feature[,i] = c(mu.i, sd.i)
          
          xx.i = (xx.i- mu.i)/sd.i
          feature.norm[,i] = xx.i
        }
        
        ## step 2: based on the constant, lag, and prediction
        ## to construct the 3D array of LSTM input
        
        ## I will temporally use 3 to predict 3 (those two tuning parameters could be adjusted)
        
        lag = 3
        prediction = 3
        
        ## for input of training data, construct x_train_arr
        
        x_train_data = list() ## this is a list of length: dim(feature)[2] + 1(for y)
        
        x_train_data[[1]] = t(sapply(
          1:(length(y.con.norm) - lag - prediction + 1), 
          function(x){
            return(y.con.norm[x:(x + lag - 1)])
          }))
        
        for (i in 1:dim(feature.norm)[2]) {
          
          xx = feature.norm[,i]
          
          x_train_data[[i+1]] = t(sapply(
            1:(length(xx) - lag - prediction + 1), 
            function(x){
              return(xx[x:(x + lag - 1)])
            }))
        }
        
        x_train_arr = array(
          data = as.numeric(unlist(x_train_data)),
          dim = c(nrow(x_train_data[[1]]),
                  lag,
                  (dim(feature.norm)[2] + 1))
        )
        
        ## for output of training data, construct y_train_arr
        
        y_train_data = t(sapply(
          (1 + lag):(length(y.con.norm) - prediction + 1),
          function(x){
            return(y.con.norm[x:(x + prediction - 1)])
          }))
        
        y_train_arr = array(
          data = as.numeric(unlist(y_train_data)),
          dim = c(
            nrow(y_train_data),
            prediction,
            1
          )
        )
        
        #################################################################
        ####################### fit lstm model ##########################
        #################################################################
        
        ## for notation simplicity, let's stick to the model output to be 'lstm_model'
        
        lstm_model <- keras::keras_model_sequential()
        lstm_model = keras::layer_lstm(lstm_model,
                                units = 50, # size of the layer
                                batch_input_shape = c(1, lag, dim(feature.norm)[2]+1), # batch size, timesteps, features
                                return_sequences = TRUE,
                                stateful = TRUE)
        lstm_model = keras::layer_dropout(lstm_model,
                                          rate = 0.5)
        lstm_model = keras::layer_lstm(lstm_model,
                                       units = 50,
                                       return_sequences = TRUE,
                                       stateful = TRUE)
        lstm_model = keras::layer_dropout(lstm_model,
                                          rate = 0.5)
        lstm_model = keras::time_distributed(lstm_model,
                                             keras::layer_dense(units = 1))
        lstm_model = keras::compile(lstm_model,
                                    loss = 'mean_squared_error',
                                    optimizer = 'adam')
        
        if(length(y.con)<=1000){
          epochs = 20
        }else{
          epochs = 50
        }
        
        keras::fit(
          x = x_train_arr,
          y = y_train_arr,
          batch_size = 1,
          epochs = epochs, 
          verbose = 0,
          shuffle = FALSE
        )
        
        #################################################################
        ##################### predict using lstm ########################
        #################################################################
        
        ## feature matrix in the post period:
        
        feature.pre = feature.norm ## a repo for feature matrix in the pre-period
        
        feature = as.matrix(feature0[post.period,]) ## this is for post period feature
        
        scale.factor.feature = matrix(NA, 2, dim(feature)[2])
        feature.norm = feature
        
        for (i in 1:dim(feature)[2]) {
          
          xx.i = feature[,i]
          mu.i = mean(xx.i)
          sd.i = stats::sd(xx.i)
          scale.factor.feature[,i] = c(mu.i, sd.i)
          
          xx.i = (xx.i- mu.i)/sd.i
          feature.norm[,i] = xx.i
        }
        
        ## the last bit of information from training data would also be used
        
        x_test_scaled = y.con.norm[(length(y.con.norm) - lag + 1):length(y.con.norm)]
        x_test_data = x_test_scaled
        
        for (i in 1:dim(feature.pre)[2]) {
          
          x_test_data = cbind(x_test_data, 
                              feature.pre[((length(y.con.norm) - lag + 1):length(y.con.norm)),i])
          
        }
        
        x_pred_arr = array(
          data = x_test_data,
          dim = c(
            1,
            lag,
            (dim(feature.pre)[2] + 1))
        )
        
        lstm_forecast = stats::predict(lstm_model,
                                       x_pred_arr, batch_size = 1)[,,1]
        
        ## Then I will need to work on the loop to move forward
        output = rep(NA, lag*ceiling(length(post.period)/lag))
        output[1:prediction] = lstm_forecast
        
        max.i = (length(output) - prediction)/lag
        
        for (i in 1:max.i) {
          
          x_test_scaled = output[(lag*(i-1) + 1):(lag*i)]
          x_test_data = x_test_scaled
          
          for (j in 1:dim(feature.pre)[2]) {
            
            x_test_data = cbind(x_test_data, 
                                feature.norm[((lag*(i-1) + 1):(lag*i)),j])
          }
          
          x_pred_arr = array(
            data = x_test_data,
            dim = c(
              1,
              lag,
              (dim(feature.pre)[2] + 1))
          )
          
          lstm_forecast = stats::predict(lstm_model,
                                         x_pred_arr, batch_size = 1)[,,1]
          
          output[(lag*i + 1):(lag*(i+1))] = lstm_forecast
          
        }
        
        output = output[1:length(post.period)]
        
        output = output*scale.factor.response[2] + scale.factor.response[1]
        
        y0 = c(y.con, output)
        
        y.con = y.con0
      }
      
      if(is.null(feature)){
        
        ## input quality control
        
        if(length(y.con) == 0){
          stop("empty vector of response in control group")
        }
        
        if(!(length(y.con) == (length(pre.period) + length(post.period)))){
          stop("length of input vector is different")
        }
        
        if(length(y.con) <= 200){
          stop("sample size is small, lstm prediction might be problematic")
        }
        
        ## make data a copy
        y.con0 = y.con
        
        ## for time series forcasting, I will only need the data from pre period
        y.con = y.con0[pre.period]
        
        ## see reference http://datasideoflife.com/?p=1171
        
        #################################################################
        ####################### data processing #########################
        #################################################################
        
        ## step 1: calculate the scale factors, for the response y.con and also feature
        ## also, normalize the response and feature
        
        scale.factor.response = c(mean(y.con), stats::sd(y.con))
        y.con.norm = (y.con - scale.factor.response[1])/scale.factor.response[2]
        
        ## step 2: based on the constant, lag, and prediction
        ## to construct the 3D array of LSTM input
        
        ## I will temporally use 3 to predict 3 (those two tuning parameters could be adjusted)
        
        lag = 3
        prediction = 3
        
        ## for input of training data, construct x_train_arr
        
        x_train_data = list() ## this is a list of length: dim(feature)[2] + 1(for y)
        
        x_train_data[[1]] = t(sapply(
          1:(length(y.con.norm) - lag - prediction + 1), 
          function(x){
            return(y.con.norm[x:(x + lag - 1)])
          }))
        
        x_train_arr = array(
          data = as.numeric(unlist(x_train_data)),
          dim = c(nrow(x_train_data[[1]]),
                  lag,
                  1)
        )
        
        ## for output of training data, construct y_train_arr
        
        y_train_data = t(sapply(
          (1 + lag):(length(y.con.norm) - prediction + 1),
          function(x){
            return(y.con.norm[x:(x + prediction - 1)])
          }))
        
        y_train_arr = array(
          data = as.numeric(unlist(y_train_data)),
          dim = c(
            nrow(y_train_data),
            prediction,
            1
          )
        )
        
        #################################################################
        ####################### fit lstm model ##########################
        #################################################################
        
        lstm_model <- keras::keras_model_sequential()
        lstm_model = keras::layer_lstm(lstm_model,
                                       units = 50, # size of the layer
                                       batch_input_shape = c(1, lag, 1), # batch size, timesteps, features
                                       return_sequences = TRUE,
                                       stateful = TRUE)
        lstm_model = keras::layer_dropout(lstm_model,
                                          rate = 0.5)
        lstm_model = keras::layer_lstm(lstm_model,
                                       units = 50,
                                       return_sequences = TRUE,
                                       stateful = TRUE)
        lstm_model = keras::layer_dropout(lstm_model,
                                          rate = 0.5)
        lstm_model = keras::time_distributed(lstm_model,
                                             keras::layer_dense(units = 1))
        lstm_model = keras::compile(lstm_model,
                                    loss = 'mean_squared_error',
                                    optimizer = 'adam')
        
        if(length(y.con)<=1000){
          epochs = 20
        }else{
          epochs = 50
        }
        
        keras::fit(lstm_model,
          x = x_train_arr,
          y = y_train_arr,
          batch_size = 1,
          epochs = epochs, 
          verbose = 0,
          shuffle = FALSE
        )
        
        
        #################################################################
        ##################### predict using lstm ########################
        #################################################################
        
        ## the last bit of information from training data would also be used
        
        x_test_scaled = y.con.norm[(length(y.con.norm) - lag + 1):length(y.con.norm)]
        x_test_data = x_test_scaled
        
        x_pred_arr = array(
          data = x_test_data,
          dim = c(
            1,
            lag,
            1)
        )
        
        lstm_forecast = stats::predict(lstm_model,
                                       x_pred_arr, batch_size = 1)[,,1]

        
        ## Then I will need to work on the loop to move forward
        output = rep(NA, lag*ceiling(length(post.period)/lag))
        output[1:prediction] = lstm_forecast
        
        max.i = (length(output) - prediction)/lag
        
        for (i in 1:max.i) {
          
          x_test_scaled = output[(lag*(i-1) + 1):(lag*i)]
          x_test_data = x_test_scaled
          
          x_pred_arr = array(
            data = x_test_data,
            dim = c(
              1,
              lag,
              1)
          )
          
          lstm_forecast = stats::predict(lstm_model,
                                         x_pred_arr, batch_size = 1)[,,1]
          
          output[(lag*i + 1):(lag*(i+1))] = lstm_forecast
          
        }
        
        output = output[1:length(post.period)]
        
        output = output*scale.factor.response[2] + scale.factor.response[1]
        
        y0 = c(y.con, output)
        
        y.con = y.con0
      }
    }
    
    ## input quality control
    
    if(length(y.exp) == 0){
      stop("empty vector of response in experiment group")
    }
    
    if(length(y.con) == 0){
      stop("empty vector of response in control group")
    }
    
    if(length(y0) == 0){
      stop("empty vector of response in control group in the absence of treatment/spillover effect")
    }
    
    if(length(pre.period) == 0){
      stop("empty time vector in pre-treatment period")
    }
    
    if(length(post.period) == 0){
      stop("empty time vector in post-treatment period")
    }
    
    if(!((length(y.exp) == length(y.con)) & (length(y.exp) == length(y0)) & (length(y.exp) == length(pre.period) + length(post.period)))){
      stop("length of input vector is different")
    }
    
    #################################################################
    ######################## estimation #############################
    #################################################################
    
    tt = c(pre.period, post.period)
    ## treatment effect model fitting
    
    dat.treatment = cbind(y.exp, y0)
    fit.treatment = CausalImpact::CausalImpact(dat.treatment, 
                                               pre.period = c(pre.period[1],pre.period[length(pre.period)]), 
                                               post.period = c(post.period[1],post.period[length(post.period)]))
    
    ## treatment effect summary table
    t.treatment = round(fit.treatment$summary[,-c(1:5,14)],3)
    x.treatment = fit.treatment$series
    
    ## spillover effect estimation
    dat.spillover = cbind(y.con, y0)
    fit.spillover = CausalImpact::CausalImpact(dat.spillover, 
                                               pre.period = c(pre.period[1],pre.period[length(pre.period)]), 
                                               post.period = c(post.period[1],post.period[length(post.period)]))
    
    ## spillover effect summary table
    t.spillover = round(fit.spillover$summary[,-c(1:5,14)],3)
    x.spillover = fit.spillover$series
    
    if(t.spillover[1,9]<=0.05){
      print(paste("Spillover effect is significant under the level of 0.05, with p-value", t.spillover[1,9]))
    }
    
    if(t.spillover[1,9]>0.05){
      print(paste("Spillover effect is NOT significant under the level of 0.05, with p-value", t.spillover[1,9]))
    }
    
    #################################################################
    ######################## graphics ###############################
    #################################################################
    
    if(is.plot == TRUE){
      
      ff = paste0(getwd(),"/SPORTSCausal_figure.pdf")
      print(paste("graphical visualization is printed as a pdf file in", ff))
      
      grDevices::pdf(file = ff, height = 16, width = 16)
      graphics::par(mfrow = c(3,3), oma = c(5,5,5,5))
      
      ## figure a: treatment response
      
      plot(tt, x.treatment$response, type = "n", 
           xlab = "time", ylab = "response of treatment group",
           main = "response of treatment group",
           ylim = range(c(as.numeric(x.treatment$response), as.numeric(x.treatment$point.pred.lower), as.numeric(x.treatment$point.pred.upper))))
      graphics::polygon(c(tt,rev(tt)),
                        c(as.numeric(x.treatment$point.pred.lower), rev(as.numeric(x.treatment$point.pred.upper))),
                        col = "gray", lty = "dashed")
      graphics::lines(x.treatment$point.pred)
      graphics::lines(x.treatment$response, col = "red", lwd = 2)
      graphics::abline(v = post.period[1], col = "blue", lty = 2, lwd = 2)
      
      ## figure b: point treatment effect
      
      plot(tt, x.treatment$point.effect,type = "n", 
           xlab = "time", ylab = "treatment effect",
           main = "treatment effect",
           ylim = range(c(as.numeric(x.treatment$point.effect.lower), as.numeric(x.treatment$point.effect.upper))))
      graphics::polygon(c(tt,rev(tt)),
                        c(as.numeric(x.treatment$point.effect.lower), rev(as.numeric(x.treatment$point.effect.upper))),
                        col = "gray", lty = "dashed")
      graphics::lines(x.treatment$point.effect, lwd = 2)
      graphics::abline(v = post.period[1], col = "blue", lty = 2, lwd = 2)
      graphics::abline(h = 0, col = "blue", lty = 2, lwd = 2)
      
      ## figure c: cumulative treatment effect
      
      plot(tt, x.treatment$cum.effect,type = "n", 
           xlab = "time", ylab = "cumulative treatment effect",
           main = "cumulative treatment effect",
           ylim = range(c(as.numeric(x.treatment$cum.effect.lower), as.numeric(x.treatment$cum.effect.upper))))
      graphics::polygon(c(tt,rev(tt)),
                        c(as.numeric(x.treatment$cum.effect.lower), rev(as.numeric(x.treatment$cum.effect.upper))),
                        col = "gray", lty = "dashed")
      graphics::lines(x.treatment$cum.effect, lwd = 2)
      graphics::abline(v = post.period[1], col = "blue", lty = 2, lwd = 2)
      graphics::abline(h = 0, col = "blue", lty = 2, lwd = 2)
      
      ## figure d: spillover response
      
      plot(tt, x.spillover$response, type = "n", 
           xlab = "time", ylab = "response of control group",
           main = "response of control group",
           ylim = range(c(as.numeric(x.spillover$response), as.numeric(x.spillover$point.pred.lower), as.numeric(x.treatment$point.pred.upper))))
      graphics::polygon(c(tt,rev(tt)),
                        c(as.numeric(x.spillover$point.pred.lower), rev(as.numeric(x.spillover$point.pred.upper))),
                        col = "gray", lty = "dashed")
      graphics::lines(x.spillover$point.pred)
      graphics::lines(x.spillover$response, col = "red", lwd = 2)
      graphics::abline(v = post.period[1], col = "blue", lty = 2, lwd = 2)
      
      ## figure e: point spillover effect
      
      plot(tt, x.spillover$point.effect,type = "n", 
           xlab = "time", ylab = "spillover effect",
           main = "spillover effect",
           ylim = range(c(as.numeric(x.spillover$point.effect.lower), as.numeric(x.spillover$point.effect.upper))))
      graphics::polygon(c(tt,rev(tt)),
                        c(as.numeric(x.spillover$point.effect.lower), rev(as.numeric(x.spillover$point.effect.upper))),
                        col = "gray", lty = "dashed")
      graphics::lines(x.spillover$point.effect, lwd = 2)
      graphics::abline(v = post.period[1], col = "blue", lty = 2, lwd = 2)
      graphics::abline(h = 0, col = "blue", lty = 2, lwd = 2)
      
      ## figure f: cumulative spillover effect
      
      plot(tt, x.spillover$cum.effect,type = "n", 
           xlab = "time", ylab = "cumulative spillover effect",
           main = "cumulative spillover effect",
           ylim = range(c(as.numeric(x.spillover$cum.effect.lower), as.numeric(x.spillover$cum.effect.upper))))
      graphics::polygon(c(tt,rev(tt)),
                        c(as.numeric(x.spillover$cum.effect.lower), rev(as.numeric(x.spillover$cum.effect.upper))),
                        col = "gray", lty = "dashed")
      graphics::lines(x.spillover$cum.effect, lwd = 2)
      graphics::abline(v = post.period[1], col = "blue", lty = 2, lwd = 2)
      graphics::abline(h = 0, col = "blue", lty = 2, lwd = 2)
      
      ## figure g: time series summary graphic
      
      plot(tt, y0, type = "n",
           xlab = "time", ylab = "response",
           main = "times series response summary",
           ylim = range(c(y0, y.exp, y.con)))
      graphics::lines(y0, lwd = 2)
      graphics::lines(y.exp, lwd = 2, col = "red")
      graphics::lines(y.con, lwd = 2, col = "blue")
      
      grDevices::dev.off()
      
    }
    
    return(list(est.treatment = t.treatment,
                est.spillover = t.spillover))
    
  }