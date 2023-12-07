# Example script to calculate vertical fluid mobility for different
# surface temperatures, EOS and depths
# Bhavik Harish Lodhia

import mobility

depths = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, \
          1.7, 1.8,1.9,2.0]


#rocks = ["Sandstone","Sandstone-clayrich","Sandstone-claypoor","Quartzite","Quartzite-quartz", \
#         "Subarkose","Subarkose-quartz","Subarkose-clayrich","Subarkose-claypoor", "Subarkose-dolomite", \
#            "Arkose","Arkose-quartzrich","Arkose-quartzpoor","Arkose-clayrich","Arkose-claypoor", "Arkose-dolomite", \
#                "Siltstone","Siltstone-orgrich","Siltstone-TOC10","Siltstone-TOC2-3","Conglomerate","Conglomerate-quartzite"] 

sandstone = ["Sandstone","Sandstone-clayrich","Sandstone-claypoor"]
quartzite = ["Quartzite","Quartzite-quartz"]
arkose = ["Subarkose","Subarkose-quartz","Subarkose-clayrich","Subarkose-claypoor", "Subarkose-dolomite", \
            "Arkose","Arkose-quartzrich","Arkose-quartzpoor","Arkose-clayrich","Arkose-claypoor", "Arkose-dolomite"]
siltstone = ["Siltstone","Siltstone-orgrich","Siltstone-TOC10","Siltstone-TOC2-3"]
conglomerate = ["Conglomerate","Conglomerate-quartzite"]

limestone = ["Limestone-OG","Limestone-WM","Micrite","Limestone-shaley","Limestone-orgrich","Limestone-TOC1-2"]
carbonates = ["Marl","Dolomite","Dolomite-sandy","Dolomite-silty","Dolomite-org"]

rocks = ["Micrite"]
# rocks = ["Micrite"]           
# equations = ["PR78", "SRK", "RK"]
equations = ["PR78"]
HCs = ["PR78", "SRK", "TWUPR", "TWUSRK"]
# tsurfs = [0., 20.]
tsurfs = [20.]

for tsurf in tsurfs:   
    for EOS in equations:
        for z in depths:
            for rock in rocks:
                out = mobility.Mobility.Mobility("H2O", rock, z, tsurf, EOS)
                print(out)      
         


# for tsurf in tsurfs:   
#     for EOS in equations:
#         for z in depths:
            # out = mobility.Mobility.Buoyancy("H2", EOS, z, tsurf)
            # print(out)