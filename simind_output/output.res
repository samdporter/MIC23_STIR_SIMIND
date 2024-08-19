


              SIMIND Monte Carlo Simulation Program    V8.0  
------------------------------------------------------------------------------
 Phantom S : h2o       Crystal...: nai       InputFile.: simind_output/outp
 Phantom B : bone      BackScatt.: pmt       OutputFile: simind_output/outp
 Collimator: pb_sb2    SourceRout: smap      SourceImg.: tmp_source        
 Cover.....: al        ScoreRout.: scattwin  DensityImg: tmp_density       
------------------------------------------------------------------------------
 PhotonEnergy.......: 208          lu177     PhotonsPerProj....: 70317          
 EnergyResolution...: 9.5          Spectra   Activity..........: 11112          
 MaxScatterOrder....: 3            g8-megp   DetectorLenght....: 19.68          
 DetectorWidth......: 25.585       SPECT     DetectorHeight....: 0.725          
 UpperEneWindowTresh: 312          BScatt    Distance to det...: 19.112         
 LowerEneWindowTresh: 104          Random    ShiftSource X.....: 0              
 PixelSize  I.......: 0.44197      Phantom   ShiftSource Y.....: 0              
 PixelSize  J.......: 0.44197      Resolut   ShiftSource Z.....: 0              
 HalfLength S.......: 28.286       Header    HalfLength P......: 28.286         
 HalfWidth  S.......: 28.286                 HalfWidth  P......: 28.286         
 HalfHeight S.......: 28.286                 HalfHeight P......: 28.286         
 SourceType.........: Integer2Map            PhantomType.......: Integer2Map  
------------------------------------------------------------------------------
 GENERAL DATA
 keV/channel........: 1                      CutoffEnergy......: 0              
 Photons/Bq.........: 0.2264                 StartingAngle.....: 360            
 CameraOffset X.....: 0                      CoverThickness....: 0              
 CameraOffset Y.....: 0                      BackscatterThickn.: 10             
 MatrixSize I.......: 128                    IntrinsicResolut..: 0.31           
 MatrixSize J.......: 128                    AcceptanceAngle...: 89.99          
 Emission type......: 3                      Initial Weight....: 35777.04646    
 NN ScalingFactor...: 0.001                  Energy Channels...: 512            
                                                                              
 SPECT DATA
 RotationMode.......: 360                    Nr of Projections.: 120            
 RotationAngle......: 3                      Projection.[start]: 1              
 Orbital fraction...: 1                      Projection...[end]: 120            
                                                                              
 COLLIMATOR DATA FOR ROUTINE: MC RayTracing       
 CollimatorCode.....: g8-megp                CollimatorType....: Parallel 
 HoleSize X.........: 0.35                   Distance X........: 0.105          
 HoleSize Y.........: 0.40415                Distance Y........: 0.29301        
 CenterShift X......: 0.2275                 X-Ray flag........: F              
 CenterShift Y......: 0.39404                CollimThickness...: 5.8            
 HoleShape..........: Hexagonal              Space Coll2Det....: 0              
 CollDepValue [57]..: 0                      CollDepValue [58].: 0              
 CollDepValue [59]..: 0                      CollDepValue [60].: 0              

 PHOTONS AFTER COLLIMATOR AND WITHIN ENER-WIN
 Geometric..........:  95.86 %          96.84 %
 Penetration........:   0.88 %           1.06 %
 Scatter in collim..:   3.27 %           2.10 %
 X-rays in collim...:   0.00 %           0.00 %
                                                                              
 IMAGE-BASED PHANTOM DATA
 RotationCentre.....:  65, 65                Bone definition...: 1170           
 CT-Pixel size......: 0.44197                Slice thickness...: 0.44197        
 StartImage.........: 1                      No of CT-Images...: 128            
 MatrixSize I.......: 128                    CTmapOrientation..: 0              
 MatrixSize J.......: 128                    StepSize..........: 0.44197        
 CenterPoint I......: 65                     ShiftPhantom X....: 0              
 CenterPoint J......: 65                     ShiftPhantom Y....: 0              
 CenterPoint K......: 65                     ShiftPhantom Z....: 0              
                                                                              
------------------------------------------------------------------------------
  Scattwin results: Window file: scattwin.win        
  
  Win  WinAdded  Range(keV)   ScaleFactor
   1       0    187.6 - 229.2    1.00
   2       1    187.6 - 229.2    1.00
  
  Win    Total    Scatter   Primary  S/P-Ratio S/T Ratio  Cps/MBq
   1   0.314E+07 0.102E+07 0.212E+07 0.482E+00 0.325E+00 0.236E+01
   2   0.314E+07 0.102E+07 0.212E+07 0.482E+00 0.325E+00 0.236E+01
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1    95.78%     1.86%     2.35%    95.92%     1.72%     2.35%
   2    95.78%     1.86%     2.35%    95.92%     1.72%     2.35%
  
  Win   SC 1  SC 2  SC 3
   1   86.4% 11.9%  1.7%
   2   86.4% 11.9%  1.7%
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 0.3031E+06     
 MaxValue projection: 1618.          
 CountRate spectrum.: 0.1905E+06     
 CountRate E-Window.: 0.9974E+05     
                                                                              
 SCATTER IN ENERGY WINDOW
 Scatter/Primary....: 1.64432        
 Scatter/Total......: 0.62183        
 Scatter order 1....: 52.73 %        
 Scatter order 2....: 32.32 %        
 Scatter order 3....: 14.95 %        
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0.4534         
 Efficiency spectrum: 0.866          
 Sensitivity Cps/MBq: 8.9758         
 Sensitivity Cpm/uCi: 19.9264        
                                                                              
 Simulation started.: 2024:08:19 15:23:19
 Simulation stopped.: 2024:08:19 15:23:40
 Elapsed time.......: 0 h, 0 m and 21 s
 DetectorHits.......: 2169           
 DetectorHits/CPUsec: 104            
                                                                              
 OTHER INFORMATION
 Compiled 2024:05:03 with INTEL Mac   
 Random number generator: Intel RAN
 Comment:EMISSION
 Energy resolution as function of 1/sqrt(E)
 Header file: simind_output/output.h00
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: simind_output\output simind_output\output /PX:0.44197402000427244/TH:0.44197402000427244/NN:0.001/CC:G8-MEGP/FI:lu177
