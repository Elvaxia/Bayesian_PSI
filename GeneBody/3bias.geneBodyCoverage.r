sample_01_3bias.Aligned.sortedByCoord.out <- c(0.0,0.0022563842300003085,0.006213629520912704,0.011764289004641933,0.018101368118685353,0.023586873650484484,0.02993309717883388,0.0382922349062869,0.045382585148785334,0.05158364110002732,0.05937925429587438,0.06913748741214218,0.0792580679452844,0.09200080928066608,0.10695878498167116,0.12072798682289898,0.13507328676540348,0.15072966710902772,0.16763426000541806,0.18494235018305974,0.20035983270294028,0.2178816735649842,0.23546866837895825,0.2526453076009515,0.26681343451627765,0.2782508107094808,0.2853514484180735,0.2957063545678064,0.3043021040154266,0.314464977464734,0.3266350498542037,0.3395298170774223,0.35313784861650727,0.3677711976096501,0.3811574771018151,0.3963474923158343,0.409011363077827,0.4214169041356757,0.4333377912353075,0.44651832140558795,0.4582563202191002,0.46767849610962325,0.47951936958407776,0.4886877879776099,0.5004497908786754,0.5092501465963919,0.5203263184245088,0.5316848240443229,0.5464862016503381,0.5597044525296307,0.5707531911148298,0.5828878288988639,0.594472658772751,0.6081412720566131,0.6224968594652118,0.6372959509676505,0.6499643939367961,0.6657659418575277,0.6774936532049458,0.6892945198668116,0.7027973906413778,0.7179062491784315,0.7329648134368024,0.7511816297861007,0.7676072839832154,0.7809912773718038,0.7967368157549114,0.8079250066582767,0.8196630054717889,0.8300224838286748,0.8409569172350492,0.8547855577692659,0.8713575225952762,0.8844043157063317,0.900555637474267,0.9138619033412547,0.9276699689432829,0.9385769691067393,0.9486289665325867,0.9609053427383635,0.968619799257245,0.9755947012691304,0.9876218921850692,0.9966885789694474,0.999893696183693,1.0,0.9885431919263966,0.9734457639072254,0.9488998698064013,0.9249872263962663,0.8959354221461712,0.8652490538388823,0.8350473394898102,0.7974477939672012,0.7506203913580712,0.6860025307166592,0.6058448810140241,0.511565969519381,0.3954319078334482,0.21458054000052582)


pdf("3bias.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,sample_01_3bias.Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()