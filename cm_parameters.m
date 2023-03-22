function [CM_stats] = cm_parameters(CM)


Accuracy = (trace(CM)/(sum(CM(:))))*100;
MissClassRate = 100-Accuracy;
Sen = diag(CM)./sum(CM,2);
Spec = (repmat(sum(CM(:)),[size(CM),1])-sum(CM,2)-sum(CM,1)+diag(CM))./(repmat(sum(CM(:)),[size(CM),1]) - sum(CM,2));

CM_stats.Accuracy               = Accuracy;
CM_stats.MissClassificationRate = MissClassRate;
CM_stats.Sensibility            = Sen;
CM_stats.Specificity            = diag(Spec);

