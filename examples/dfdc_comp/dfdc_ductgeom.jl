

duct_coords = [
    0.304542 0.159526
    0.296440 0.162842
    0.285240 0.167270
    0.273301 0.171796
    0.261206 0.176188
    0.249026 0.180416
    0.236794 0.184470
    0.224549 0.188335
    0.212269 0.192017
    0.200020 0.195490
    0.187807 0.198746
    0.175647 0.201771
    0.163586 0.204547
    0.151644 0.207051
    0.139859 0.209267
    0.128380 0.211153
    0.117508 0.212618
    0.106950 0.213649
    0.096610 0.214255
    0.086499 0.214421
    0.076647 0.214136
    0.067113 0.213390
    0.057935 0.212174
    0.049193 0.210487
    0.040949 0.208332
    0.033312 0.205736
    0.026403 0.202735
    0.020339 0.199393
    0.015227 0.195790
    0.011135 0.192004
    0.008090 0.188086
    0.006112 0.184067
    0.005242 0.180005
    0.005484 0.176154
    0.006854 0.172546
    0.009324 0.169289
    0.012842 0.166404
    0.017419 0.163862
    0.023110 0.161648
    0.029956 0.159771
    0.037937 0.158256
    0.046983 0.157103
    0.057025 0.156294
    0.067995 0.155792
    0.079836 0.155546
    0.092531 0.155498
    0.106044 0.155585
    0.120000 0.155721
    0.134222 0.155902
    0.148679 0.156177
    0.163490 0.156523
    0.178507 0.156897
    0.193399 0.157258
    0.208123 0.157586
    0.222751 0.157864
    0.237332 0.158088
    0.251898 0.158254
    0.266505 0.158365
    0.281130 0.158423
    0.294972 0.158441
    0.304466 0.158439
]

duct_coords = reverse(duct_coords, dims=1)

hub_coords = [
    0.306379 0.035928
    0.296701 0.039168
    0.285768 0.041795
    0.274184 0.043576
    0.262095 0.044564
    0.248578 0.044839
    0.232723 0.044847
    0.216580 0.044866
    0.200407 0.044882
    0.184232 0.044898
    0.168060 0.044913
    0.151899 0.044930
    0.135773 0.044950
    0.120000 0.044952
    0.109361 0.044999
    0.099680 0.044922
    0.090039 0.044561
    0.080630 0.043916
    0.071451 0.042966
    0.062540 0.041690
    0.053933 0.040075
    0.045705 0.038107
    0.037900 0.035771
    0.030604 0.033068
    0.023901 0.030001
    0.017880 0.026585
    0.012632 0.022848
    0.008231 0.018825
    0.004736 0.014551
    0.002179 0.010047
    0.000586 0.005293
    0.000000 0.000000
]

hub_coords = reverse(hub_coords, dims=1)
hub_coords[end,1] = maximum(duct_coords[:,1])