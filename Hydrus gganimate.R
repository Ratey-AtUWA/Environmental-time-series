library(reshape2)
library(gganimate)
library(gifski)
library(av)
library(stringr)

bentonite <- read.csv("bentonite.csv")

solute0$Step <- as.factor(paste0(as.character(solute0$Time)," days"))
solute0$Step <- factor(solute0$Step, levels = paste0(seq(0,60,6)," days"))
with(solute0, plot(Conc ~ Depth, type="l", xlim=c(0,-300), ylim=c(0,140), 
                   subset = solute0$Step == levels(solute0$Step)[1]))
# plot(solute0$Conc ~ solute0$Depth, type="l", xlim=c(0,-300), ylim=c(0,140), 
#      col = rainbow(11)[solute0$Step], xaxt="n", xlab="")

for(i in 2:nlevels(solute0$Step)) {
  with(solute0, lines(Conc ~ Depth, 
                      col = rainbow(nlevels(solute0$Step),v=0.8,end=0.9)[i], lwd=2, 
                      subset = solute0$Step==levels(solute0$Step)[i]))
  }

# _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._
#
# make the animation
require(gganimate)

solute0$Dist <- -1*solute0$Depth

bentAnim <- ggplot(solute0, aes(x=Dist, y=Conc, col=Step)) + 
  geom_line(size=0.75) +
  scale_color_manual(values=rainbow(nlevels(solute0$Step), v=0.8, end=0.9)) +
  theme_bw() +
  theme(legend.position="none",
        title = element_text(size=14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold", colour="black"),
        axis.text = element_text(size = 14),
        panel.border = element_rect(colour = 1,fill=NA)) + 
  labs(y="Concentration", x = "Distance (cm)") + 
  transition_states(solute0$Step,
                    transition_length = 2,
                    state_length = 1,
                    wrap = FALSE) +
  ease_aes() +
  enter_fade() + 
  exit_fade() + 
  ease_aes("cubic-in") +
  ggtitle("You are now looking at approximately {closest_state}")

#  shadow_trail(distance = 0.2) +
  
require(gifski)

anim_save("bentAnim2.gif", animate(bentAnim, renderer = gifski_renderer(), 
                                  duration = 20))

require(av)
options(gganimate.dev_args = list(width = 800, height = 600))
anim_save("bentAnim.mp4", animate(bentAnim, renderer = av_renderer(), 
                                  duration = 20))


# geom_text(aes(x = 5, y = 0.35, label = "Smectite\n+ Mg"), size = 4, color="red2") +
#   geom_text(aes(x=7, y=0.65, label="Smectite + Mg + Gly"), size=4, color="blue3") +