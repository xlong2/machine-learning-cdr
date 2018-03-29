
aab=ggplot(all_frames_I_care, aes(mean_simu_error_count, error_count,colour = factor(loop_type)))+
  geom_point( position=position_dodge(width = 0.90))+
  scale_colour_manual(breaks = all_frames_I_care$loop_type, name="loop type",
                      values =col_vec)+
  facet_wrap(~sig_stat)+xlab(" error count from random assignment")+ylab("blindBLAST error count")+
  ggtitle("Misclassifications categorized by the significance test")+   
  theme_classic()+theme_bw()+
  theme(strip.text.x = element_text(size = 11),
        axis.text.x = element_text(size = 11, colour = "black"),plot.title =element_text(hjust = 0.5,size=18,face="bold"),
        plot.margin=unit(c(2,2,2,2),"mm"),
        axis.text.y = element_text(size=11)) 
ggsave(aab,units = c("in"),width=9,height=8,file="./proline_classifier/Plots/misclassifications_categorized_by_significance.pdf")
ggsave(aab,units = c("in"),dpi=300,width=9,height=8,file="./proline_classifier/Plots/misclassifications_categorized_by_significance.png")
