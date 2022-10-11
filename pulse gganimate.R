library(reshape2)
library(gganimate)
library(gifski)
library(av)
library(stringr)

npts <- 500
nd <- 11
ds <- data.frame(d=seq(1,96,l=nd), s=seq(0.3,6,l=nd))
x0 <- runif(npts)
y0 <- rnorm(npts,ds$d[1],ds$s[1])
plot(x0,y0,ylim=c(100,0),pch=16, xaxt="n", xlab="",
     ylab = "Depth", cex = 2, 
     col = rainbow(nd, v=0.7, end=0.8, alpha=0.25)[1])
for(i in 2:nd) {
  points(x0, rnorm(npts, ds$d[i], ds$s[i]), cex = 2,
         pch = 16, col = rainbow(nd, v=0.7, end=0.8, alpha=0.25)[i])
}
pulse <- data.frame(x0 = runif(npts))
for(i in 0:(nd-1)) {
  pulse[,paste0("time_",i)] <- rnorm(npts, ds$d[i+1], ds$s[i+1])
}
str(pulse)

pulseFact <- melt(pulse, measure.vars=2:ncol(pulse), 
                     variable.name = "Timestep", value.name = "Particles")
str(pulseFact)


####  ------------ Optionally repeat last timestep at end  ---------------- ####
end <- pulseFact[nrow(pulseFact)-(npts-1):nrow(pulseFact),]
end$Timestep <- as.character(end$Timestep)
end$Timestep <- rep("time_end",npts)
pulseFact$Timestep <- as.character(pulseFact$Timestep)
pulseFact <- rbind(pulseFact,end)
pulseFact$Timestep <- as.factor(pulseFact$Timestep)
levels(pulseFact$Timestep)
pulseFact$Timestep <- factor(pulseFact$Timestep,
                             levels=c(paste0("time_",seq(0,10,1)),"time_end"))



# _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._
#
# make the animation ####
require(gganimate)

pulseAnim <- ggplot(pulseFact, aes(x=x0, y=Particles)) + 
  scale_y_reverse() +
  geom_point(size=4, shape = 16, color="#0000B040") +
  theme_bw() +
  theme(legend.position="none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        title = element_text(size=18, face = "bold", colour = "blue3"),
        axis.title = element_text(size = 22, face = "bold", colour="black"),
        axis.text = element_text(size = 18),
        panel.border = element_rect(colour = 1,fill=NA)) + 
  labs(y="Depth (cm)", x = "") + 
  transition_states(pulseFact$Timestep,
                    transition_length = 1,
                    state_length = c(3,rep(0,10),3),
                    wrap = FALSE) +
  enter_fade() + 
  exit_fade() + 
  shadow_wake(wake_length = 0.05, size = T, colour="#60606040") +
  ggtitle("Distribution of solute at {closest_state}")

#  shadow_trail(distance = 0.2) +
  
require(gifski)
options(gganimate.dev_args = list(width = 600, height = 800))
anim_save("pulseAnim.gif", animate(pulseAnim, renderer = gifski_renderer(), 
                                  duration = 20))

require(av)
options(gganimate.dev_args = list(width = 400, height = 800))
anim_save("bentAnim.mp4", animate(bentAnim, renderer = av_renderer(), 
                                  duration = 20))


# geom_text(aes(x = 5, y = 0.35, label = "Smectite\n+ Mg"), size = 4, color="red2") +
#   geom_text(aes(x=7, y=0.65, label="Smectite + Mg + Gly"), size=4, color="blue3") +