sample_01_nobias.Aligned.sortedByCoord.out <- c(0.0,0.10898175984520528,0.20891430169476793,0.2976491684885406,0.37865185736663787,0.4464954829493384,0.503555252382293,0.5619774119123158,0.6128998123643042,0.6565340871513954,0.6927630820627455,0.723374343605742,0.7539553718040047,0.7811219216313913,0.8058112267746501,0.8189929650785972,0.830304015177139,0.840349043964952,0.8529128882859016,0.8651233803902747,0.8746412152292916,0.884512402284885,0.894585774833386,0.9077183839521407,0.9207206117717307,0.927579801858217,0.931608395043999,0.9401190815865703,0.9430668326981182,0.9511561419984619,0.9594022867746123,0.9650313576472411,0.9683097859668152,0.9739048443266184,0.9760646388910409,0.9822889287381169,0.9857147446132682,0.990276200499984,0.9876647953486,0.9924341054803606,0.9941970873951517,0.9943803770476005,1.0,0.9967688112815726,0.9960998985293368,0.9944030520561509,0.995597269173137,0.9927269910074695,0.9936302121813925,0.9920542990871419,0.9891821313374287,0.9887305207504672,0.9833093041228834,0.9810380240997549,0.9799363966010162,0.9786722648743332,0.9761232159964627,0.9763613035862415,0.9696230468786906,0.9628016484731217,0.9530136031155462,0.9449110667268814,0.9382937433982658,0.9332258789872585,0.928413108422443,0.9210682952361696,0.9145964698790855,0.9053545143107647,0.8981911011928944,0.8870520032425262,0.8805934049737631,0.8755047551382514,0.87038398237396,0.8598760054949104,0.850879695852552,0.8414204381189568,0.8348843669043133,0.8226058497742892,0.8115328872655262,0.8043902595721604,0.7917621693936514,0.778958347898877,0.767469676900024,0.7566253540608106,0.7424440257966014,0.7265639614751604,0.7046750198878721,0.6826802615940153,0.6552019303990613,0.6265898487765889,0.5967344208519378,0.5598459611085812,0.5193200520769363,0.4762942233526134,0.427884969681624,0.37281871141705575,0.3071840095839703,0.23293658366983677,0.14720993467707955,0.049316254013004115)


pdf("nobias.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,sample_01_nobias.Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()