
# head(cross.table)
#     Species Clustered.with Num..reference.genes Num..predicted.clusters C.... S.... D.... F.... M.... Rand.Score Adj..Rand.Score       NMI        AMI
# 1 C. maxima            All                26988                   22676 98.49 84.09 14.40  1.03  0.47  0.9999950       0.9516716 0.9963985        0.9534170
# 5 C. maxima           Self                26988                   26584 89.77 77.17 12.60  3.48  6.75  0.9999979       0.9767026 0.9984210        0.9776109
# 6 C. maxima         Others                26988                   25042 89.59 80.22  9.37  3.61  6.79  0.9999944       0.9395084 0.9956583        0.9407056
# 7 C. maxima  C. reticulata                26988                   25094 89.64 79.71  9.93  3.57  6.79  0.9999935       0.9295787 0.9949561        0.9320180
# 8 C. maxima    C. sinensis                26988                   25116 89.69 79.97  9.72  3.53  6.79  0.9999941       0.9365153 0.9954255        0.9381042
# 9 C. maxima     C. hindsii                26988                   25111 89.59 79.79  9.80  3.61  6.79  0.9999934       0.9289454 0.9949648        0.9319704

# This function call will produce a barplot with ylim cutoff that doesn't mess with adobe illustrator
# The key is the combination of coord_cartesian(ylim=<...>) with scale_y_continuous(limits=<...>, oob = scales::squish)
# Most other details are not relevant for bringing about the adobe friendly behaviour

ggplot(cross.table, aes(x = Clustered.with, y = NMI, fill = Clustered.with)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw() +
  ggtitle("Normalised Mutual Information") +
  coord_cartesian(ylim=c(0.89, 1)) +
  scale_y_continuous(limits=c(0.89, 1), expand = expansion(mult = c(0, 0.05)), oob = scales::squish) +
  facet_wrap(Species ~ ., scales='free') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_text(size=12),
        strip.text.x = element_text(size=13),
        legend.text = element_text(size=12),
        legend.title = element_text(size=13))