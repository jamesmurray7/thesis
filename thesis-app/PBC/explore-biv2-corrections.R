rm(list=ls())
source(".Rprofile")
library(ggplot2)
load(save.dir.file.path(x = "Joint/Multivs/biv2.RData"))
PBCredx$histologic2 <- ifelse(PBCredx$histologic%in%c("3","4"), 1, 0)

# Correction: QQplot for residuals ----------------------------------------
ggplot(out.all, aes(sample = pr)) + 
  geom_qq_line(lty = 5, col = 'red', lwd = 0.25) + 
  geom_qq(size=.005, alpha = .2) + 
  facet_wrap(~response, scales = 'free') +
  labs(y = 'Sample quantiles', x = 'Theoretical quantiles') + 
  theme_csda()+
  theme(
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5),
    axis.title = element_text(size=5),
    strip.text = element_text(size=7,vjust=1)
  )

ggsave("biv2-qqnorm-corrections.png", width = 140, height = 60, units = "mm")


# Correction: QQplot for REs ----------------------------------------------
# Posterior distns
a <- cond.ranefs(joint.biv2)
a$qnames <- c(
  expression(b[10]), expression(b[11]), 
  expression(b[12]), expression(b[13]),
  expression(b[20]), expression(b[21])
)
png('biv2-postb-corrections.png', width = 140, height = 130, units = 'mm',
    res = 1000L)
par(mar = c(2, 4, 2, 3))
plot(a)
dev.off()

# QQnorms
REs <- ranef(joint.biv2)
colnames(REs) <- c('b[10]', 'b[11]', 'b[12]', 'b[13]',
                   'b[20]', 'b[21]')
REs.long <- tidyr::pivot_longer(as.data.frame(REs), everything())
REs.long$name <- forcats::fct_inorder(REs.long$name)

ggplot(REs.long, aes(sample = value)) + 
  geom_qq_line(lty = 5, col = 'red', lwd = 0.25) + 
  geom_qq(size=.005, alpha = .2) + 
  facet_wrap(~name, scales = 'free', labeller = label_parsed, nrow = 3, ncol = 2) +
  labs(y = 'Sample quantiles', x = 'Theoretical quantiles') + 
  theme_csda()+
  theme(
    axis.text.x = element_text(size=5),
    axis.text.y = element_text(size=5),
    axis.title = element_text(size=5),
    strip.text = element_text(size=7,vjust=1)
  )
ggsave("biv2-qqnorm-REs-corrections.png", width = 140, height = 90, units = "mm")

