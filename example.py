# Example script to calculate vertical fluid mobility for different
# surface temperatures, EOS and depths
# Bhavik Harish Lodhia
# UNSW Sydney

import mobility

depths = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, \
          1.7, 1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2, \
          3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]

rocks = ["Sandstone"] 
# rocks = ["Micrite"]           
equations = ["PR78", "SRK", "RK"]
HCs = ["PR78", "SRK", "TWUPR", "TWUSRK"]
tsurfs = [0., 20.]

for tsurf in tsurfs:   
    for EOS in equations:
        for z in depths:
            for rock in rocks:
                out = mobility.Mobility.Mobility("CO2", rock, z, tsurf, EOS)
                print(out)
                # out = mobility.Mobility.Mobility("methane", rock, z, tsurf, EOS)
                # print(out)
                # out = mobility.Mobility.Mobility("H2", rock, z, tsurf, EOS)
                # print(out)      
                # out = mobility.Mobility.Mobility("drygas", rock, z, tsurf, EOS)
                # print(out)  
                # out = mobility.Mobility.Mobility("wetgas", rock, z, tsurf, EOS)
                # print(out)                
                

# for tsurf in tsurfs:   
#     for EOS in HCs:
#         for z in depths:
#             for rock in rocks:
                # out = mobility.Mobility.Mobility("mixture", rock, z, tsurf, EOS)
                # print(out)                    
                # out = mobility.Mobility.Mobility("voil", rock, z, tsurf, EOS)
                # print(out)                
                # out = mobility.Mobility.Mobility("lightoil", rock, z, tsurf, EOS)
                # print(out) 
                # out = mobility.Mobility.Mobility("blackoil", rock, z, tsurf, EOS)
                # print(out)  


# for tsurf in tsurfs:   
#     for EOS in equations:
#         for z in depths:
#             out = mobility.Mobility.Buoyancy("CO2", EOS, z, tsurf)
#             print(out)
            # out = mobility.Mobility.Buoyancy("methane", EOS, z, tsurf)
            # print(out)   
            # out = mobility.Mobility.Buoyancy("H2", EOS, z, tsurf)
            # print(out)    
            # out = mobility.Mobility.Buoyancy("mixture", EOS, z, tsurf)
            # print(out)
            # out = mobility.Mobility.Buoyancy("drygas", EOS, z, tsurf)
            # print(out)   
            # out = mobility.Mobility.Buoyancy("wetgas", EOS, z, tsurf)
            # print(out)    
            # out = mobility.Mobility.Buoyancy("voil", EOS, z, tsurf)
            # print(out)   
            # out = mobility.Mobility.Buoyancy("lightoil", EOS, z, tsurf)
            # print(out)    
            # out = mobility.Mobility.Buoyancy("blackoil", EOS, z, tsurf)
            # print(out)   