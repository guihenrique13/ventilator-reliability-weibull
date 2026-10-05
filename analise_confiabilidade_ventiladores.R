# ============================================================
#  ANÁLISE DE CONFIABILIDADE — VENTILADORES PULMONARES
#  Macrorregião de Saúde do Vale do Jequitinhonha
#  Padrão editorial: tons de azul, Times New Roman, 600 DPI
# ============================================================
#
#  Pacotes necessários:
#    install.packages(c("fitdistrplus","ggplot2","dplyr","tidyr",
#                       "patchwork","scales","boot","nortest"))
#
#  R version: 4.4.1
#  Autores: Caique Menezes de Abreu, Leandra Cardoso Moreira,
#           Guilherme Henrique Fernandes Carneiro, Paulo Cesar de Resende Andrade
#  Data: 2025
#  Dados: idênticos a data/tbf_ventiladores.csv (vetor embutido abaixo)
# ============================================================

library(fitdistrplus)   # v. 1.2-1  — ajuste MLE das distribuições
library(ggplot2)        # v. 3.5.1  — geração de gráficos
library(dplyr)          # v. 1.3.2  — manipulação de dados
library(tidyr)          # v. 1.3.1  — reorganização de dados
library(patchwork)      # v. 1.3.0  — composição de painéis
library(scales)         # v. 1.3.0  — formatação de eixos
library(boot)           # v. 1.3-30 — intervalos de confiança bootstrap
library(nortest)        # v. 1.0-4  — teste de Lilliefors

set.seed(2024)   # reprodutibilidade

# ============================================================
# 1. DADOS — Tempos entre falhas (TBF, dias)
# ============================================================
tbf <- c(
  13.95, 48.17, 27.76, 21.74,  7.08,  7.08,  3.53, 36.81, 21.84, 26.54,
   1.75, 53.30, 34.02,  8.89,  7.92,  7.97, 11.75, 18.97, 15.80, 11.35,
  22.27,  6.53, 11.38, 13.69, 16.60, 30.78,  8.49, 18.59, 21.50,  3.03,
  22.10,  7.55,  3.82, 47.78, 51.94, 32.29, 11.76,  5.07, 25.40, 16.07,
   5.93, 17.93,  2.47, 41.42, 10.34, 24.41, 11.98, 18.80, 19.76,  8.02,
  53.19, 30.17, 45.95, 39.70, 21.71, 43.12,  4.73,  8.37,  2.98, 12.41,
  14.40, 10.73, 33.74, 13.39, 11.03, 19.62,  6.58, 31.88,  4.20, 61.42,
  30.00,  8.46,  0.72, 32.78, 26.48, 27.60, 29.94,  4.18, 13.44,  5.72,
  36.54, 22.74, 12.58,  3.76, 11.96, 12.40, 27.63, 23.33, 38.87, 17.14,
   5.85, 26.80, 29.33, 20.30, 29.92, 17.88, 18.90, 15.66,  2.01,  5.44
)
n <- length(tbf)
cat("N =", n, "| min =", min(tbf), "| max =", max(tbf), "\n")

# ============================================================
# 2. AJUSTE MLE DAS DISTRIBUIÇÕES (loc/threshold = 0 fixo)
# ============================================================
fit_weibull  <- fitdist(tbf, "weibull",  method = "mle")
fit_lognorm  <- fitdist(tbf, "lnorm",    method = "mle")
fit_expon    <- fitdist(tbf, "exp",      method = "mle")
fit_normal   <- fitdist(tbf, "norm",     method = "mle")

# Parâmetros Weibull
beta_w <- fit_weibull$estimate["shape"]
eta_w  <- fit_weibull$estimate["scale"]
cat(sprintf("Weibull: beta = %.4f | eta = %.4f dias\n", beta_w, eta_w))

# MTTF Weibull
mttf_w <- eta_w * gamma(1 + 1/beta_w)
cat(sprintf("MTTF Weibull = %.4f dias\n", mttf_w))

# Percentis
t10_w <- qweibull(0.10, shape = beta_w, scale = eta_w)
t50_w <- qweibull(0.50, shape = beta_w, scale = eta_w)
t90_w <- qweibull(0.90, shape = beta_w, scale = eta_w)
cat(sprintf("T10=%.4f | T50=%.4f | T90=%.4f\n", t10_w, t50_w, t90_w))

# Intervalos de confiança bootstrap (B = 5000)
stat_beta <- function(data, i) {
  fit <- fitdist(data[i], "weibull", method = "mle")
  fit$estimate["shape"]
}
stat_eta <- function(data, i) {
  fit <- fitdist(data[i], "weibull", method = "mle")
  fit$estimate["scale"]
}
boot_beta <- boot(tbf, stat_beta, R = 5000)
boot_eta  <- boot(tbf, stat_eta,  R = 5000)
ci_beta   <- boot.ci(boot_beta, type = "perc")$percent[4:5]
ci_eta    <- boot.ci(boot_eta,  type = "perc")$percent[4:5]
cat(sprintf("IC 95%% beta: [%.4f; %.4f]\n", ci_beta[1], ci_beta[2]))
cat(sprintf("IC 95%% eta:  [%.4f; %.4f]\n", ci_eta[1],  ci_eta[2]))

# Critérios de informação
ll_w  <- logLik(fit_weibull)[1]
ll_ln <- logLik(fit_lognorm)[1]
ll_e  <- logLik(fit_expon)[1]
ll_n  <- logLik(fit_normal)[1]
aic_w  <- -2*ll_w  + 2*2;  bic_w  <- -2*ll_w  + 2*log(n)
aic_ln <- -2*ll_ln + 2*2;  bic_ln <- -2*ll_ln + 2*log(n)
aic_e  <- -2*ll_e  + 2*1;  bic_e  <- -2*ll_e  + 1*log(n)
aic_n  <- -2*ll_n  + 2*2;  bic_n  <- -2*ll_n  + 2*log(n)

# ============================================================
# 3. PALETA E TEMA PADRÃO (azul monocromático, Times New Roman)
# ============================================================
AZ1 <- "#1A3A5C"   # azul escuro  — Weibull / principal
AZ2 <- "#2B6CB0"   # azul médio   — Lognormal
AZ3 <- "#63B3ED"   # azul claro   — Exponencial
AZ4 <- "#A0AEC0"   # cinza-azul   — Normal
AZ_HIST <- "#BEE3F8"  # azul pálido — histograma

COR_DIST <- c("Weibull" = AZ1, "Lognormal" = AZ2,
              "Exponencial" = AZ3, "Normal" = AZ4)
LT_DIST  <- c("Weibull" = "solid", "Lognormal" = "dashed",
              "Exponencial" = "dotted", "Normal" = "dotdash")

tema <- theme_bw() +
  theme(
    text             = element_text(family = "serif", size = 11),
    plot.title       = element_text(family = "serif", size = 11, hjust = 0.5),
    axis.title       = element_text(family = "serif", size = 10),
    axis.text        = element_text(family = "serif", size = 9),
    legend.text      = element_text(family = "serif", size = 9),
    legend.title     = element_text(family = "serif", size = 9, face = "bold"),
    panel.grid.major = element_line(color = "#cccccc", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "#333333", linewidth = 0.8),
    legend.background = element_rect(fill = "white", color = "#cccccc",
                                      linewidth = 0.4),
    plot.margin      = margin(8, 8, 8, 8)
  )

# Função auxiliar para salvar
salvar <- function(p, nome, w = 7, h = 4.5) {
  ggsave(filename = paste0("figs/", nome, ".png"),
         plot = p, width = w, height = h,
         dpi = 600, bg = "white")
  ggsave(filename = paste0("figs/", nome, ".pdf"),
         plot = p, width = w, height = h, bg = "white")
  message("Salvo: figs/", nome, ".png / .pdf")
}

dir.create("figs", showWarnings = FALSE)

# ============================================================
# 4. GRADE DE PONTOS PARA AS CURVAS
# ============================================================
t_seq <- seq(0.1, 65, by = 0.2)

curvas <- data.frame(
  t = t_seq,
  # PDFs
  f_weibull  = dweibull(t_seq, shape = beta_w, scale = eta_w),
  f_lognorm  = dlnorm  (t_seq, meanlog = fit_lognorm$estimate["meanlog"],
                                sdlog   = fit_lognorm$estimate["sdlog"]),
  f_expon    = dexp    (t_seq, rate = fit_expon$estimate["rate"]),
  f_normal   = dnorm   (t_seq, mean = fit_normal$estimate["mean"],
                                sd   = fit_normal$estimate["sd"]),
  # CDFs
  F_weibull  = pweibull(t_seq, shape = beta_w, scale = eta_w),
  F_lognorm  = plnorm  (t_seq, meanlog = fit_lognorm$estimate["meanlog"],
                                sdlog   = fit_lognorm$estimate["sdlog"]),
  F_expon    = pexp    (t_seq, rate = fit_expon$estimate["rate"]),
  F_normal   = pnorm   (t_seq, mean = fit_normal$estimate["mean"],
                                sd   = fit_normal$estimate["sd"]),
  # Confiabilidade
  R_weibull  = 1 - pweibull(t_seq, shape = beta_w, scale = eta_w),
  R_lognorm  = 1 - plnorm  (t_seq, meanlog = fit_lognorm$estimate["meanlog"],
                                     sdlog   = fit_lognorm$estimate["sdlog"]),
  R_expon    = 1 - pexp    (t_seq, rate = fit_expon$estimate["rate"]),
  # Hazard
  h_weibull  = (beta_w / eta_w) * (t_seq / eta_w)^(beta_w - 1),
  h_lognorm  = {
    mu <- fit_lognorm$estimate["meanlog"]
    sg <- fit_lognorm$estimate["sdlog"]
    f  <- dlnorm(t_seq, mu, sg)
    S  <- 1 - plnorm(t_seq, mu, sg)
    ifelse(S > 1e-10, f / S, 0)
  },
  # Risco acumulado
  H_weibull  = (t_seq / eta_w)^beta_w,
  H_lognorm  = -log(pmax(1 - plnorm(t_seq,
                                     fit_lognorm$estimate["meanlog"],
                                     fit_lognorm$estimate["sdlog"]), 1e-15))
)

# ============================================================
# FIGURA 1 — Histograma + curvas de densidade
# ============================================================
df_hist <- data.frame(x = tbf)

df_dens <- tidyr::pivot_longer(
  curvas[, c("t","f_weibull","f_lognorm","f_expon","f_normal")],
  cols = -t, names_to = "dist", values_to = "densidade"
) |>
  mutate(dist = recode(dist,
    f_weibull = "Weibull", f_lognorm = "Lognormal",
    f_expon   = "Exponencial", f_normal = "Normal"))

fig1 <- ggplot() +
  geom_histogram(data = df_hist, aes(x = x, y = after_stat(density)),
                 bins = 12, fill = AZ_HIST, color = AZ1,
                 linewidth = 0.5, alpha = 0.85) +
  geom_line(data = df_dens,
            aes(x = t, y = densidade, color = dist, linetype = dist),
            linewidth = 0.9) +
  scale_color_manual(values = COR_DIST, name = "Distribuição") +
  scale_linetype_manual(values = LT_DIST, name = "Distribuição") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Figura 1a — Função de densidade de probabilidade f(t)",
    x = "Tempo entre falhas (dias)",
    y = "Densidade de probabilidade f(t)"
  ) +
  tema +
  theme(legend.position = "right")

salvar(fig1, "fig1_histograma_densidade", w = 7, h = 4.5)

# ============================================================
# FIGURA 2 — Curva de confiabilidade R(t)
# ============================================================
df_Rt <- tidyr::pivot_longer(
  curvas[, c("t","R_weibull","R_lognorm","R_expon")],
  cols = -t, names_to = "dist", values_to = "Rt"
) |>
  mutate(dist = recode(dist,
    R_weibull = "Weibull", R_lognorm = "Lognormal", R_expon = "Exponencial"))

fig2 <- ggplot() +
  geom_ribbon(data = filter(df_Rt, dist == "Weibull"),
              aes(x = t, ymin = 0, ymax = Rt), fill = AZ1, alpha = 0.10) +
  geom_line(data = df_Rt,
            aes(x = t, y = Rt, color = dist, linetype = dist), linewidth = 1.0) +
  geom_vline(xintercept = c(t10_w, t50_w, t90_w),
             linetype = "dashed", color = "#888888", linewidth = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             color = "#888888", linewidth = 0.6) +
  annotate("text", x = t10_w + 0.8, y = 0.04, label = "T(10)",
           size = 3, color = "#666666", hjust = 0) +
  annotate("text", x = t50_w + 0.8, y = 0.04, label = "T(50)",
           size = 3, color = "#666666", hjust = 0) +
  annotate("text", x = t90_w + 0.8, y = 0.04, label = "T(90)",
           size = 3, color = "#666666", hjust = 0) +
  scale_color_manual(values = c("Weibull"=AZ1,"Lognormal"=AZ2,"Exponencial"=AZ3),
                     name = "Modelo") +
  scale_linetype_manual(values = c("Weibull"="solid","Lognormal"="dashed",
                                   "Exponencial"="dotted"), name = "Modelo") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Figura 1b — Função de confiabilidade R(t)",
    x = "Tempo (dias)",
    y = "Confiabilidade R(t)"
  ) +
  tema

salvar(fig2, "fig2_confiabilidade_Rt", w = 7, h = 4.5)

# ============================================================
# FIGURA 3 — Função de taxa de falha h(t)
# ============================================================
df_ht <- tidyr::pivot_longer(
  curvas[, c("t","h_weibull","h_lognorm")],
  cols = -t, names_to = "dist", values_to = "ht"
) |>
  mutate(dist = recode(dist,
    h_weibull = "Weibull — monotonicamente crescente",
    h_lognorm = "Lognormal — não monótona"))

fig3 <- ggplot() +
  geom_ribbon(data = filter(df_ht, dist == "Weibull — monotonicamente crescente"),
              aes(x = t, ymin = 0, ymax = ht), fill = AZ1, alpha = 0.10) +
  geom_line(data = df_ht,
            aes(x = t, y = ht, color = dist, linetype = dist), linewidth = 1.0) +
  scale_color_manual(
    values = c("Weibull — monotonicamente crescente" = AZ1,
               "Lognormal — não monótona"            = AZ2),
    name = "Modelo") +
  scale_linetype_manual(
    values = c("Weibull — monotonicamente crescente" = "solid",
               "Lognormal — não monótona"            = "dashed"),
    name = "Modelo") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, NA)) +
  labs(
    title = "Figura 1c — Função de taxa de falha h(t) (hazard function)",
    x = "Tempo (dias)",
    y = "Taxa de falha h(t)  [falhas/dia]"
  ) +
  tema +
  theme(legend.position = "bottom")

salvar(fig3, "fig3_hazard_ht", w = 7, h = 4.5)

# ============================================================
# FIGURA 4 — Painel Figura 1 completo (f(t) + R(t) + h(t))
# ============================================================
fig_painel <- fig1 + fig2 + fig3 +
  plot_layout(ncol = 3) &
  theme(
    legend.position = "bottom",
    plot.title      = element_text(size = 9)
  )

salvar(fig_painel, "fig4_painel_fRh", w = 14, h = 4.5)

# ============================================================
# FIGURA 5 — Ajuste das distribuições: histograma + CDF
#            (equivalente à Figura 2 do artigo)
# ============================================================
df_Ft <- tidyr::pivot_longer(
  curvas[, c("t","F_weibull","F_lognorm","F_expon","F_normal")],
  cols = -t, names_to = "dist", values_to = "Ft"
) |>
  mutate(dist = recode(dist,
    F_weibull = "Weibull", F_lognorm = "Lognormal",
    F_expon   = "Exponencial", F_normal = "Normal"))

# ECDF empírica
df_ecdf <- data.frame(
  x = sort(tbf),
  F_emp = (seq_along(tbf) - 0.375) / (n + 0.25)
)

p5a <- ggplot() +
  geom_histogram(data = df_hist, aes(x = x, y = after_stat(density)),
                 bins = 12, fill = AZ_HIST, color = AZ1,
                 linewidth = 0.5, alpha = 0.85) +
  geom_line(data = df_dens,
            aes(x = t, y = densidade, color = dist, linetype = dist),
            linewidth = 0.8) +
  scale_color_manual(values = COR_DIST, name = NULL) +
  scale_linetype_manual(values = LT_DIST, name = NULL) +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "TBF (dias)", y = "Densidade") +
  tema + theme(legend.position = "right", plot.title = element_blank())

p5b <- ggplot() +
  geom_step(data = df_ecdf, aes(x = x, y = F_emp, color = "ECDF empírica"),
            linewidth = 0.7) +
  geom_line(data = df_Ft,
            aes(x = t, y = Ft, color = dist, linetype = dist),
            linewidth = 0.8) +
  scale_color_manual(
    values = c("ECDF empírica" = "#333333", COR_DIST),
    name = NULL) +
  scale_linetype_manual(
    values = c("ECDF empírica" = "solid", LT_DIST),
    name = NULL) +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(x = "TBF (dias)", y = "F(t) — Probabilidade acumulada") +
  tema + theme(legend.position = "right", plot.title = element_blank())

fig5 <- p5a + p5b +
  plot_annotation(
    title = "Figura 2 — Ajuste de distribuições probabilísticas aos tempos de falha",
    theme = theme(plot.title = element_text(family = "serif", size = 11, hjust = 0.5))
  )

salvar(fig5, "fig5_ajuste_distribuicoes", w = 12, h = 4.8)

# ============================================================
# FIGURA 6 — Comparação R(t), h(t), H(t): Weibull vs. Lognormal
#            (Figura 4 do artigo)
# ============================================================
df_Ht <- tidyr::pivot_longer(
  curvas[, c("t","H_weibull","H_lognorm")],
  cols = -t, names_to = "dist", values_to = "Ht"
) |>
  mutate(dist = recode(dist,
    H_weibull = "Weibull", H_lognorm = "Lognormal"))

p6a <- ggplot() +
  geom_line(data = filter(df_Rt, dist != "Exponencial"),
            aes(x = t, y = Rt, color = dist, linetype = dist), linewidth = 1.0) +
  scale_color_manual(values = c("Weibull"=AZ1,"Lognormal"=AZ2), name = "Modelo") +
  scale_linetype_manual(values = c("Weibull"="solid","Lognormal"="dashed"), name = "Modelo") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 20)) +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(title = "Confiabilidade R(t)", x = "Tempo (dias)", y = "R(t)") +
  tema + theme(legend.position = "bottom")

p6b <- ggplot() +
  geom_line(data = df_ht,
            aes(x = t, y = ht, color = dist, linetype = dist), linewidth = 1.0) +
  scale_color_manual(
    values = c("Weibull — monotonicamente crescente" = AZ1,
               "Lognormal — não monótona" = AZ2),
    labels = c("Weibull","Lognormal"), name = "Modelo") +
  scale_linetype_manual(
    values = c("Weibull — monotonicamente crescente" = "solid",
               "Lognormal — não monótona" = "dashed"),
    labels = c("Weibull","Lognormal"), name = "Modelo") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Taxa de falha h(t)", x = "Tempo (dias)", y = "h(t)  [falhas/dia]") +
  tema + theme(legend.position = "bottom")

p6c <- ggplot() +
  geom_line(data = df_Ht,
            aes(x = t, y = Ht, color = dist, linetype = dist), linewidth = 1.0) +
  scale_color_manual(values = c("Weibull"=AZ1,"Lognormal"=AZ2), name = "Modelo") +
  scale_linetype_manual(values = c("Weibull"="solid","Lognormal"="dashed"), name = "Modelo") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Risco acumulado H(t)", x = "Tempo (dias)", y = "H(t)") +
  tema + theme(legend.position = "bottom")

fig6 <- p6a + p6b + p6c +
  plot_annotation(
    title = "Figura 4 — Comparação dos modelos Weibull e Lognormal: R(t), h(t) e H(t)",
    theme = theme(plot.title = element_text(family = "serif", size = 11, hjust = 0.5))
  ) &
  theme(legend.position = "bottom")

salvar(fig6, "fig6_comparacao_modelos_RhH", w = 13, h = 4.5)

# ============================================================
# FIGURA 7 — QQ-Plots para as quatro distribuições
# ============================================================
qq_weibull <- function() {
  p_i <- (seq_along(tbf) - 0.375) / (n + 0.25)
  q_emp <- sort(tbf)
  q_teo <- qweibull(p_i, shape = beta_w, scale = eta_w)
  data.frame(teo = q_teo, emp = q_emp, dist = "Weibull")
}
qq_lognorm <- function() {
  p_i <- (seq_along(tbf) - 0.375) / (n + 0.25)
  q_emp <- sort(tbf)
  q_teo <- qlnorm(p_i, fit_lognorm$estimate["meanlog"],
                       fit_lognorm$estimate["sdlog"])
  data.frame(teo = q_teo, emp = q_emp, dist = "Lognormal")
}
qq_expon <- function() {
  p_i <- (seq_along(tbf) - 0.375) / (n + 0.25)
  q_emp <- sort(tbf)
  q_teo <- qexp(p_i, rate = fit_expon$estimate["rate"])
  data.frame(teo = q_teo, emp = q_emp, dist = "Exponencial")
}
qq_norm <- function() {
  p_i <- (seq_along(tbf) - 0.375) / (n + 0.25)
  q_emp <- sort(tbf)
  q_teo <- qnorm(p_i, fit_normal$estimate["mean"],
                       fit_normal$estimate["sd"])
  data.frame(teo = q_teo, emp = q_emp, dist = "Normal")
}

df_qq <- bind_rows(qq_weibull(), qq_lognorm(), qq_expon(), qq_norm()) |>
  mutate(dist = factor(dist, levels = c("Weibull","Lognormal","Exponencial","Normal")))

fig7 <- ggplot(df_qq, aes(x = teo, y = emp, color = dist)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#333333", linewidth = 0.8) +
  facet_wrap(~dist, scales = "free", ncol = 2) +
  scale_color_manual(values = COR_DIST, guide = "none") +
  labs(
    title = "Figura — Gráficos quantil-quantil (QQ-plots) para as quatro distribuições",
    x = "Quantis teóricos",
    y = "Quantis empíricos"
  ) +
  tema

salvar(fig7, "fig7_qqplots", w = 9, h = 7)

# ============================================================
# FIGURA 8 — PP-Plots para as quatro distribuições
# ============================================================
pp_data <- function(dist_name, cdf_fn) {
  p_emp <- (seq_along(tbf) - 0.375) / (n + 0.25)
  p_teo <- cdf_fn(sort(tbf))
  data.frame(teo = p_teo, emp = p_emp, dist = dist_name)
}

df_pp <- bind_rows(
  pp_data("Weibull",
          function(x) pweibull(x, shape=beta_w, scale=eta_w)),
  pp_data("Lognormal",
          function(x) plnorm(x, fit_lognorm$estimate["meanlog"],
                                 fit_lognorm$estimate["sdlog"])),
  pp_data("Exponencial",
          function(x) pexp(x, rate = fit_expon$estimate["rate"])),
  pp_data("Normal",
          function(x) pnorm(x, fit_normal$estimate["mean"],
                                fit_normal$estimate["sd"]))
) |>
  mutate(dist = factor(dist, levels = c("Weibull","Lognormal","Exponencial","Normal")))

fig8 <- ggplot(df_pp, aes(x = teo, y = emp, color = dist)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#333333", linewidth = 0.8) +
  facet_wrap(~dist, ncol = 2) +
  scale_color_manual(values = COR_DIST, guide = "none") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Figura — Gráficos PP-Plot para as quatro distribuições",
    x = "Probabilidades teóricas",
    y = "Probabilidades empíricas"
  ) +
  tema

salvar(fig8, "fig8_ppplots", w = 9, h = 7)

# ============================================================
# FIGURA 9 — Gráfico de probabilidade Weibull (linearizado)
# ============================================================
p_i   <- (seq_along(tbf) - 0.375) / (n + 0.25)
x_wei <- log(sort(tbf))
y_wei <- log(-log(1 - p_i))
lm_wei <- lm(y_wei ~ x_wei)
R2 <- summary(lm_wei)$r.squared

df_wei_prob <- data.frame(x = x_wei, y = y_wei)
df_wei_line <- data.frame(
  x = range(x_wei),
  y = predict(lm_wei, newdata = data.frame(x_wei = range(x_wei)))
)

fig9 <- ggplot(df_wei_prob, aes(x = x, y = y)) +
  geom_point(color = AZ1, size = 2, alpha = 0.75) +
  geom_line(data = df_wei_line, aes(x = x, y = y),
            color = AZ2, linewidth = 1.2) +
  annotate("label",
           x = min(x_wei) + 0.3 * diff(range(x_wei)),
           y = max(y_wei) - 0.15 * diff(range(y_wei)),
           label = sprintf("β̂ = %.3f\nη̂ = %.2f dias\nR² = %.4f",
                           beta_w, eta_w, R2),
           family = "serif", size = 3.5,
           fill = "white", color = "#333333", label.size = 0.3) +
  labs(
    title = "Figura — Gráfico de probabilidade Weibull",
    x = "ln(t)",
    y = "ln(−ln(1−F))"
  ) +
  tema

salvar(fig9, "fig9_weibull_probplot", w = 6.5, h = 5)

# ============================================================
# FIGURA 10 — Curva da Banheira interpretada
# ============================================================
df_bath <- data.frame(t = t_seq, h = (beta_w/eta_w)*(t_seq/eta_w)^(beta_w-1))

fig10 <- ggplot() +
  annotate("rect", xmin=0,  xmax=6,  ymin=0, ymax=Inf,
           fill="#BEE3F8", alpha=0.45) +
  annotate("rect", xmin=6,  xmax=17, ymin=0, ymax=Inf,
           fill="#C6F6D5", alpha=0.40) +
  annotate("rect", xmin=17, xmax=65, ymin=0, ymax=Inf,
           fill="#FED7D7", alpha=0.35) +
  annotate("text", x=3.0, y=0.24, label="Mortalidade\ninfantil\n(β<1)",
           family="serif", size=3.2, color="#2A4365", hjust=0.5, vjust=1) +
  annotate("text", x=11.5,y=0.24, label="Vida útil\n(β≈1)",
           family="serif", size=3.2, color="#22543D", hjust=0.5, vjust=1) +
  annotate("text", x=41,   y=0.24, label="Desgaste / Envelhecimento  (β>1)",
           family="serif", size=3.2, color="#742A2A", hjust=0.5, vjust=1) +
  geom_line(data = df_bath, aes(x = t, y = h),
            color = AZ1, linewidth = 2.2) +
  geom_vline(xintercept = c(t10_w, t50_w, t90_w),
             linetype = "dashed", color = "#718096", linewidth = 0.8) +
  annotate("text", x=t10_w+0.5, y=0.005, label="T(10)", size=3,
           color="#718096", hjust=0, family="serif") +
  annotate("text", x=t50_w+0.5, y=0.005, label="T(50)", size=3,
           color="#718096", hjust=0, family="serif") +
  annotate("text", x=t90_w+0.5, y=0.005, label="T(90)", size=3,
           color="#718096", hjust=0, family="serif") +
  scale_x_continuous(limits=c(0,65), breaks=seq(0,60,10)) +
  scale_y_continuous(limits=c(0,0.27),
                     expand=expansion(mult=c(0,0))) +
  labs(
    title = "Figura — Curva da Banheira: posicionamento dos ventiladores no ciclo de vida",
    x = "Tempo de operação sem intervenção (dias)",
    y = "Taxa de falha h(t)  [falhas/dia]"
  ) +
  tema

salvar(fig10, "fig10_curva_banheira", w = 8, h = 4.8)

# ============================================================
# FIGURA 11 — Bootstrap: distribuição de β e η
# ============================================================
df_boot <- data.frame(
  beta_b = boot_beta$t,
  eta_b  = boot_eta$t
)

p11a <- ggplot(df_boot, aes(x = beta_b)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = AZ_HIST, color = AZ1, linewidth = 0.4, alpha = 0.85) +
  geom_vline(xintercept = beta_w,
             color = AZ1, linewidth = 1.5, linetype = "solid") +
  geom_vline(xintercept = ci_beta,
             color = AZ2, linewidth = 1.0, linetype = "dashed") +
  annotate("label", x = beta_w, y = Inf, vjust = 1.3,
           label = sprintf("MLE = %.4f", beta_w),
           family = "serif", size = 3.2, fill="white", color=AZ1) +
  annotate("label",
           x = mean(range(df_boot$beta_b)), y = Inf, vjust = 3.0,
           label = sprintf("IC 95%% [%.3f; %.3f]", ci_beta[1], ci_beta[2]),
           family = "serif", size = 3.2, fill="white", color=AZ2) +
  labs(title = "β (forma)", x = "β", y = "Densidade") +
  tema

p11b <- ggplot(df_boot, aes(x = eta_b)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = AZ_HIST, color = AZ1, linewidth = 0.4, alpha = 0.85) +
  geom_vline(xintercept = eta_w,
             color = AZ1, linewidth = 1.5, linetype = "solid") +
  geom_vline(xintercept = ci_eta,
             color = AZ2, linewidth = 1.0, linetype = "dashed") +
  annotate("label", x = eta_w, y = Inf, vjust = 1.3,
           label = sprintf("MLE = %.3f dias", eta_w),
           family = "serif", size = 3.2, fill="white", color=AZ1) +
  annotate("label",
           x = mean(range(df_boot$eta_b)), y = Inf, vjust = 3.0,
           label = sprintf("IC 95%% [%.2f; %.2f]", ci_eta[1], ci_eta[2]),
           family = "serif", size = 3.2, fill="white", color=AZ2) +
  labs(title = "η (escala, dias)", x = "η (dias)", y = "Densidade") +
  tema

fig11 <- p11a + p11b +
  plot_annotation(
    title = sprintf(
      "Figura — Bootstrap (B=5.000): distribuição de β̂ e η̂\nIC 95%% β: [%.4f; %.4f]  |  IC 95%% η: [%.2f; %.2f] dias",
      ci_beta[1], ci_beta[2], ci_eta[1], ci_eta[2]
    ),
    theme = theme(plot.title = element_text(family="serif", size=10, hjust=0.5))
  )

salvar(fig11, "fig11_bootstrap_ci", w = 11, h = 4.5)

# ============================================================
# FIGURA 12 — Curva g(T) vs T (modelo de substituição por idade)
# ============================================================
Cf <- 6693.00
Cp <- 1249.10

gT <- function(T_val) {
  Rt_val <- 1 - pweibull(T_val, shape = beta_w, scale = eta_w)
  Ft_val <- pweibull(T_val, shape = beta_w, scale = eta_w)
  Mt_val <- integrate(function(u)
    1 - pweibull(u, shape = beta_w, scale = eta_w),
    lower = 0, upper = T_val)$value
  (Cp * Rt_val + Cf * Ft_val) / Mt_val
}

T_range <- seq(1, 40, by = 0.2)
gT_vals <- sapply(T_range, gT)
# Ótimo contínuo (a grade de 0,2 dia arredondava T* para 14,6)
opt_gT  <- optimize(gT, interval = c(1, 40))
T_star  <- opt_gT$minimum
gT_star <- opt_gT$objective
gT_rtf  <- Cf / mttf_w  # custo assintótico RTF

df_gT <- data.frame(T = T_range, gT = gT_vals)

fig12 <- ggplot(df_gT, aes(x = T, y = gT)) +
  geom_ribbon(aes(ymin = gT_rtf, ymax = gT),
              fill = AZ_HIST, alpha = 0.4) +
  geom_line(color = AZ1, linewidth = 1.6) +
  geom_hline(yintercept = gT_rtf, linetype = "longdash",
             color = "#718096", linewidth = 0.9) +
  geom_vline(xintercept = T_star, linetype = "dashed",
             color = AZ1, linewidth = 1.0) +
  geom_vline(xintercept = 12, linetype = "dotted",
             color = AZ2, linewidth = 0.9) +
  geom_vline(xintercept = 14, linetype = "dotdash",
             color = AZ3, linewidth = 0.9) +
  annotate("label", x = T_star, y = gT_star,
           label = sprintf("T* = %.1f dias\ng(T*) = R$ %.0f/dia",
                           T_star, gT_star),
           vjust = -0.3, family = "serif", size = 3.2,
           fill = "white", color = AZ1, label.size = 0.3) +
  annotate("text", x = 37, y = gT_rtf + 1.5,
           label = sprintf("RTF: g(∞) = R$ %.0f/dia", gT_rtf),
           family = "serif", size = 3.0, color = "#718096", hjust = 1) +
  annotate("text", x = 12, y = gT_vals[which(T_range == 12)] + 3,
           label = "T=12\n(Prev.)",
           family = "serif", size = 2.8, color = AZ2, hjust = 0.5) +
  annotate("text", x = 14, y = gT_vals[which.min(abs(T_range-14))] + 3,
           label = "T=14\n(PMBC)",
           family = "serif", size = 2.8, color = AZ3, hjust = 0.5) +
  scale_x_continuous(limits = c(1, 40), breaks = seq(0, 40, 5)) +
  scale_y_continuous(
    labels = scales::label_dollar(prefix="R$ ", big.mark=".", decimal.mark=",")) +
  labs(
    title = sprintf(
      "Figura 5 — Função de custo de ciclo g(T) vs. intervalo de manutenção T\nCf = R$ %.0f  |  Cp = R$ %.0f  |  Weibull (β=%.4f, η=%.3f dias)",
      Cf, Cp, beta_w, eta_w),
    x = "Intervalo de manutenção preventiva T (dias)",
    y = "Custo médio por unidade de tempo g(T)  [R$/dia]"
  ) +
  tema

salvar(fig12, "fig12_gT_custo_ciclo", w = 8, h = 4.8)

# ============================================================
# FIGURA 13 — Custos anuais por estratégia (barras comparativas)
# ============================================================
df_custos <- data.frame(
  estrategia = c("Corretiva\n(RTF)", "Preventiva\n(tempo)",
                 "PMBC", "RCM", "Preditiva"),
  custo = c(125109.88, 113396.36, 112646.73, 102838.00, 76467.70),
  tipo  = c("RTF", "Prev", "PMBC", "RCM", "Pred")
)
df_custos$estrategia <- factor(df_custos$estrategia,
                                levels = df_custos$estrategia)

cores_barras <- c("RTF"="#1A3A5C","Prev"="#2B6CB0","PMBC"="#2C7BB6",
                  "RCM"="#63B3ED","Pred"="#BEE3F8")

fig13 <- ggplot(df_custos, aes(x = estrategia, y = custo, fill = tipo)) +
  geom_col(color = "black", linewidth = 0.5, width = 0.65) +
  geom_text(aes(label = scales::dollar(custo, prefix="R$ ", big.mark=".",
                                       decimal.mark=",")),
            vjust = -0.5, family = "serif", size = 3.2) +
  scale_fill_manual(values = cores_barras, guide = "none") +
  scale_y_continuous(
    limits = c(0, 145000),
    labels = scales::label_dollar(prefix="R$ ", big.mark=".", decimal.mark=",")) +
  labs(
    title = "Figura — Custo anual por estratégia de manutenção (R$/equipamento/ano)",
    x = NULL,
    y = "Custo anual estimado (R$/equipamento/ano)"
  ) +
  tema

salvar(fig13, "fig13_custos_anuais_barras", w = 8, h = 5)

# ============================================================
# FIGURA 14 — Economia potencial (cenário atual vs. proposto)
# ============================================================
df_eco <- data.frame(
  estrategia = rep(c("PMBC","RCM","Preditiva"), each = 2),
  cenario    = rep(c("Atual (RTF)","Proposto"), 3),
  custo      = c(125110, 112647,
                 125110, 102838,
                 125110,  76468),
  cor        = c(AZ1, AZ2, AZ1, AZ3, AZ1, AZ_HIST)
)
df_eco$estrategia <- factor(df_eco$estrategia, levels=c("PMBC","RCM","Preditiva"))
df_eco$cenario    <- factor(df_eco$cenario,
                             levels=c("Atual (RTF)","Proposto"))

df_anotacao <- data.frame(
  estrategia = c("PMBC","RCM","Preditiva"),
  economia   = c(125110-112647, 125110-102838, 125110-76468),
  pct        = c((125110-112647)/125110*100,
                 (125110-102838)/125110*100,
                 (125110-76468)/125110*100),
  y_pos      = c(112647, 102838, 76468)
)

fig14 <- ggplot(df_eco,
                aes(x = estrategia, y = custo,
                    fill = interaction(cenario, estrategia))) +
  geom_col(position = position_dodge(width = 0.72),
           width = 0.68, color = "black", linewidth = 0.5) +
  geom_text(
    data = df_eco |> filter(cenario == "Atual (RTF)"),
    aes(label = scales::dollar(custo, prefix="R$ ", big.mark=".", decimal.mark=",")),
    position = position_dodge(width = 0.72),
    vjust = -0.4, hjust = 0.5, family = "serif", size = 3.0) +
  geom_text(
    data = df_eco |> filter(cenario == "Proposto"),
    aes(label = scales::dollar(custo, prefix="R$ ", big.mark=".", decimal.mark=",")),
    position = position_dodge(width = 0.72),
    vjust = 1.4, hjust = 0.5, family = "serif", size = 3.0,
    color = "white", fontface = "bold") +
  geom_label(
    data = df_anotacao,
    aes(x = estrategia, y = y_pos + 8000,
        label = sprintf("−R$ %s\n(−%.1f%%)",
                        formatC(economia, format="f", digits=0,
                                big.mark=".", decimal.mark=","), pct)),
    inherit.aes = FALSE,
    family = "serif", size = 3.0, color = AZ1,
    fill = "white", label.size = 0.3, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "Atual (RTF).PMBC"      = AZ1, "Proposto.PMBC"     = AZ2,
      "Atual (RTF).RCM"       = AZ1, "Proposto.RCM"      = AZ3,
      "Atual (RTF).Preditiva" = AZ1, "Proposto.Preditiva"= "#63B3ED"
    ),
    breaks = c("Atual (RTF).PMBC","Proposto.PMBC"),
    labels = c("Cenário atual (Corretiva/RTF)","Cenário proposto"),
    name = NULL) +
  scale_y_continuous(
    limits = c(0, 148000),
    labels = scales::label_dollar(prefix="R$ ", big.mark=".", decimal.mark=",")) +
  labs(
    title = paste0("Figura — Economia potencial: cenário atual vs. cenário proposto\n",
                   "(R$/equipamento/ano)"),
    x = NULL,
    y = "Custo anual (R$/equipamento/ano)"
  ) +
  tema +
  theme(legend.position = "top")

salvar(fig14, "fig14_economia_potencial", w = 8.5, h = 5.5)

# ============================================================
# FIGURA 15 — Custo acumulado 5 anos
# ============================================================
anos <- 0:5
df_acum <- data.frame(
  ano = rep(anos, 5),
  custo = c(
    125109.88 * anos,
    1100  + 113396.36 * anos,
    1500  + 112646.73 * anos,
    3800  + 102838.00 * anos,
    7200  +  76467.70 * anos
  ),
  estrategia = rep(
    c("Corretiva (RTF)","Preventiva (tempo)","PMBC","RCM","Preditiva"), each=6)
)
df_acum$estrategia <- factor(df_acum$estrategia,
  levels=c("Corretiva (RTF)","Preventiva (tempo)","PMBC","RCM","Preditiva"))

cores_acum <- c("Corretiva (RTF)"=AZ1,"Preventiva (tempo)"=AZ2,
                "PMBC"="#2C7BB6","RCM"=AZ3,"Preditiva"=AZ_HIST)
formas_acum <- c("Corretiva (RTF)"=16,"Preventiva (tempo)"=15,
                 "PMBC"=17,"RCM"=18,"Preditiva"=25)

fig15 <- ggplot(df_acum, aes(x = ano, y = custo,
                              color = estrategia,
                              shape = estrategia,
                              group = estrategia)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = cores_acum, name = "Estratégia") +
  scale_shape_manual(values = formas_acum, name = "Estratégia") +
  scale_x_continuous(breaks = 0:5) +
  scale_y_continuous(
    labels = scales::label_dollar(prefix="R$ ", big.mark=".", decimal.mark=",")) +
  labs(
    title = "Figura — Custo acumulado em horizonte de 5 anos por estratégia de manutenção",
    x = "Tempo (anos)",
    y = "Custo acumulado (R$/equipamento)"
  ) +
  tema +
  theme(legend.position = "right")

salvar(fig15, "fig15_custo_acumulado_5anos", w = 9, h = 5)

# ============================================================
# FIGURA 16 — F(t) + custo operacional esperado h(t)*Cf
# ============================================================
df_ftcusto <- data.frame(
  t     = t_seq,
  Ft    = pweibull(t_seq, shape=beta_w, scale=eta_w),
  custo = (beta_w/eta_w)*(t_seq/eta_w)^(beta_w-1) * Cf
)

fig16 <- ggplot(df_ftcusto, aes(x = t)) +
  geom_line(aes(y = Ft, color = "F(t)"), linewidth = 1.2) +
  geom_line(aes(y = custo / max(custo), color = "h(t)·Cf"),
            linewidth = 1.2, linetype = "dashed") +
  scale_color_manual(
    values = c("F(t)"=AZ1, "h(t)·Cf"=AZ2),
    labels = c("F(t) — Probabilidade acumulada de falha",
               "Custo operacional esperado h(t)·Cf (escala normalizada)"),
    name = NULL) +
  scale_x_continuous(limits=c(0,65), breaks=seq(0,60,10)) +
  scale_y_continuous(
    name = "F(t) — Probabilidade acumulada de falha",
    limits = c(0, 1.05),
    sec.axis = sec_axis(
      ~ . * max(df_ftcusto$custo),
      name = "Custo operacional esperado (R$/dia)",
      labels = scales::label_dollar(prefix="R$ ", big.mark=".", decimal.mark=","))) +
  labs(
    title = "Figura — F(t) e custo operacional esperado h(t)·Cf em função do tempo",
    x = "Tempo de operação sem intervenção (dias)"
  ) +
  tema +
  theme(legend.position = "bottom")

salvar(fig16, "fig16_Ft_custo_operacional", w = 8, h = 5)

# ============================================================
# TABELA DE RESULTADOS (console)
# ============================================================
cat("\n", strrep("=", 65), "\n")
cat("TABELA 3 — RESULTADOS COMPLETOS (para revisão)\n")
cat(strrep("=", 65), "\n")
cat(sprintf("%-14s %8s %8s %8s %8s %7s %7s\n",
            "Distribuição","LogLik","AIC","BIC","T10","T50","T90"))
cat(strrep("-", 65), "\n")
resultados <- list(
  list("Weibull",   ll_w,  aic_w,  bic_w,
       qweibull(0.1,beta_w,eta_w), qweibull(0.5,beta_w,eta_w), qweibull(0.9,beta_w,eta_w)),
  list("Lognormal", ll_ln, aic_ln, bic_ln,
       qlnorm(0.1,fit_lognorm$estimate[1],fit_lognorm$estimate[2]),
       qlnorm(0.5,fit_lognorm$estimate[1],fit_lognorm$estimate[2]),
       qlnorm(0.9,fit_lognorm$estimate[1],fit_lognorm$estimate[2])),
  list("Exponencial",ll_e, aic_e,  bic_e,
       qexp(0.1,fit_expon$estimate[1]), qexp(0.5,fit_expon$estimate[1]),
       qexp(0.9,fit_expon$estimate[1])),
  list("Normal",    ll_n,  aic_n,  bic_n,
       qnorm(0.1,fit_normal$estimate[1],fit_normal$estimate[2]),
       qnorm(0.5,fit_normal$estimate[1],fit_normal$estimate[2]),
       qnorm(0.9,fit_normal$estimate[1],fit_normal$estimate[2]))
)
for (r in resultados) {
  cat(sprintf("%-14s %8.2f %8.2f %8.2f %8.3f %8.3f %8.3f\n",
              r[[1]], r[[2]], r[[3]], r[[4]], r[[5]], r[[6]], r[[7]]))
}
cat("\nWeibull:\n")
cat(sprintf("  β̂ = %.6f  IC95%% [%.4f; %.4f]\n", beta_w, ci_beta[1], ci_beta[2]))
cat(sprintf("  η̂ = %.6f  IC95%% [%.4f; %.4f] dias\n", eta_w, ci_eta[1], ci_eta[2]))
cat(sprintf("  MTTF = %.4f dias\n", mttf_w))
cat(sprintf("  T* ótimo g(T) = %.1f dias  |  g(T*) = R$ %.2f/dia\n",
            T_star, gT_star))
cat(sprintf("  g(RTF) = R$ %.2f/dia\n", gT_rtf))
cat("\nFiguras geradas em: ./figs/\n")
