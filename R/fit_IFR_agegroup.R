df <- data.table::fread(here::here("ifr_data.csv"))

df[, agemid := (agehigh + agelow) / 2] 

fit <- lm(data = df, boot::logit(middle / 100) ~ agemid)

preds <- predict(fit, newdata = data.frame(agemid = seq(2, 95, 1)), se.fit = TRUE)
preds <- data.table(x = seq(2, 95, 1), 
                    y = preds$fit, 
                    ymin = preds$fit - 1.96 * preds$se.fit,
                    ymax = preds$fit + 1.96 * preds$se.fit)

p1 <- df %>%
  ggplot(aes(x = agemid, 
             y = log(middle / 100), 
             col = author,
             ymin = log(bottom / 100),
             ymax = log(top / 100))) + 
  geom_point() +
  geom_line() +
  geom_errorbar() +
  geom_line(inherit.aes = FALSE, data = preds, aes(x = x, y = y), col = "blue", size = 1.2) +
  geom_ribbon(inherit.aes = FALSE, data = preds, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.2) +
  cowplot::theme_cowplot() +
  labs(y = "logit(Infection Fatality Ratio)",
       x = "Age group mid-point")

p2 <- df %>%
  ggplot(aes(x = agemid,
             y = middle,
             col = author,
             ymin = bottom,
             ymax = top)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  geom_line(inherit.aes = FALSE, data = preds, 
            aes(x = x, y = boot::inv.logit(y)* 100), col = "blue", size = 1.2) +
  geom_ribbon(inherit.aes = FALSE, data = preds, 
              aes(x = x, ymin = boot::inv.logit(ymin) * 100, 
                  ymax = boot::inv.logit(ymax) * 100), alpha = 0.2) +
  cowplot::theme_cowplot() +
  labs(y = "Infection Fatality Ratio (%)",
       x = "Age group mid-point")

p1 + p2





