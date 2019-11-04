class task1:    
    float_tlu_time = [5014.0000000, 18076.0000000, 49121.0000000, 107285.0000000, 207577.0000000, 337898.0000000, 678805.0000000, 1102976.0000000, 1754699.0000000, 2589886.0000000, 3527384.0000000, 4857919.0000000, 6809750.0000000, 8931249.0000000, 11754777.0000000, 14921192.0000000, 18700242.0000000, 22971065.0000000, 28392056.0000000, 34854213.0000000]
    float_optimized_tlu_time = [3021.0000000, 11029.0000000, 27050.0000000, 55147.0000000, 100048.0000000, 169451.0000000, 309824.0000000, 501505.0000000, 746454.0000000, 1124999.0000000, 1553098.0000000, 2100624.0000000, 2868628.0000000, 3591022.0000000, 4481424.0000000, 5639026.0000000, 6833681.0000000, 8254989.0000000, 9907317.0000000, 12302720.0000000]
    float_triangular_time = [9728.0000000, 23038.0000000, 39978.0000000, 63327.0000000, 93736.0000000, 130909.0000000, 176065.0000000, 234766.0000000, 304484.0000000, 433831.0000000, 661648.0000000, 1008418.0000000, 1457165.0000000, 2006738.0000000, 2624403.0000000, 3353131.0000000, 4227934.0000000, 5207866.0000000, 6313109.0000000, 7583256.0000000]
    
    double_tlu_time = [3801.0000000, 18048.0000000, 48159.0000000, 108345.0000000, 203256.0000000, 357977.0000000, 656746.0000000, 1078894.0000000, 1725617.0000000, 2569833.0000000, 3649700.0000000, 4970754.0000000, 6881321.0000000, 9096223.0000000, 11977891.0000000, 15239543.0000000, 18666159.0000000, 22755563.0000000, 28840215.0000000, 34671138.0000000]
    double_optimized_tlu_time = [2006.0000000, 11028.0000000, 28067.0000000, 56149.0000000, 101696.0000000, 171249.0000000, 299797.0000000, 461735.0000000, 741984.0000000, 1127024.0000000, 1541104.0000000, 2164397.0000000, 2802435.0000000, 3585504.0000000, 4451806.0000000, 5453532.0000000, 6688814.0000000, 8059440.0000000, 9722858.0000000, 12314783.0000000]
    double_triangular_time = [10414.0000000, 22058.0000000, 46175.0000000, 64193.0000000, 92537.0000000, 134318.0000000, 177492.0000000, 238270.0000000, 301543.0000000, 461265.0000000, 639657.0000000, 1002439.0000000, 1423462.0000000, 1992788.0000000, 2622463.0000000, 3383357.0000000, 4201136.0000000, 5209086.0000000, 6354864.0000000, 7543308.0000000]
    
    long_double_tlu_time = [4036.0000000, 17076.0000000, 52171.0000000, 105280.0000000, 201557.0000000, 344948.0000000, 628698.0000000, 1120006.0000000, 1749680.0000000, 2594814.0000000, 3653686.0000000, 5032411.0000000, 6840225.0000000, 8880675.0000000, 11688092.0000000, 14662602.0000000, 18065554.0000000, 22832175.0000000, 27631536.0000000, 34443089.0000000]
    long_double_optimized_tlu_time = [2393.0000000, 10997.0000000, 27065.0000000, 56149.0000000, 101269.0000000, 172461.0000000, 298799.0000000, 501305.0000000, 751006.0000000, 1137992.0000000, 1551911.0000000, 2145707.0000000, 2789417.0000000, 3535311.0000000, 4433817.0000000, 5434434.0000000, 6637649.0000000, 8004256.0000000, 9709807.0000000, 12153321.0000000]
    long_double_triangular_time = [10344.0000000, 23061.0000000, 41126.0000000, 63362.0000000, 106158.0000000, 144025.0000000, 170384.0000000, 247013.0000000, 347889.0000000, 418075.0000000, 795810.0000000, 1165991.0000000, 1410583.0000000, 1970702.0000000, 2728345.0000000, 3341013.0000000, 4205082.0000000, 5177033.0000000, 6311445.0000000, 7471729.0000000]
    
    single_thread_optimized_tlu_time = [117312.0000000, 980608.0000000, 5572816.0000000, 17824955.0000000, 38911648.0000000, 73314981.0000000, 127444424.0000000, 197835440.0000000, 287596839.0000000, 404761569.0000000]
    multi_thread_optimized_tlu_time = [40106.0000000, 325542.0000000, 1300432.0000000, 3391023.0000000, 6680734.0000000, 12075664.0000000, 20721405.0000000, 32036427.0000000, 48503511.0000000, 70161006.0000000]
    
    single_thread_triangular_time = [30080.0000000, 258688.0000000, 1583241.0000000, 4498966.0000000, 9123794.0000000, 15837120.0000000, 24523716.0000000, 36551544.0000000, 51911257.0000000, 69949165.0000000]
    
class task2:
    float_fault = [1.2031860352, 2.0000000000, 2.0000000000, 2.0000000000, 3.0000000000, 3.0000000000, 3.0000000000, 3.0000000000, 3.0000000000, 3.4999847412, 3.5012817383, 4.0000000000, 4.0000000000, 4.1250000000, 4.2500000000, 5.0000000000, 5.1874694824, 6.0000000000, 6.0000000000, 6.0000000000, 6.5000000000, 6.8398437500, 6.9375000000, 6.9531250000, 7.0000000000, 7.0000000000, 7.0000000000, 7.5000000000, 7.6562500000, 7.9061889648, 8.0000000000, 8.0000000000, 8.0000000000, 8.0000000000, 8.7500000000, 9.0000000000, 10.0000000000, 10.0001220703, 10.0050048828, 10.0625000000, 11.0000000000, 11.6250000000, 12.0000000000, 12.0000000000, 12.0000000000, 12.4998779297, 13.0000000000, 13.0000000000, 13.0000610352, 13.1718750000, 13.5000000000, 14.0000305176, 14.7500000000, 14.9997558594, 15.0000000000, 15.0000000000, 15.0000000000, 15.4997558594, 15.5000000000, 15.5000000000, 15.9999389648, 16.0000000000, 16.3749694824, 16.9995117188, 17.0000000000, 17.0000305176, 17.2500000000, 17.5000000000, 18.0000000000, 18.0000000000, 18.0000000000, 19.7500000000, 21.0000000000, 21.0000000000, 21.9999542236, 22.0000000000, 22.0000000000, 22.1250000000, 22.7500000000, 22.9998779297, 22.9999618530, 23.0000000000, 24.0000000000, 24.0000000000, 24.0000000000, 25.0000000000, 25.0000000000, 26.0000000000, 26.0000000000, 26.3750000000, 27.0000000000, 27.5000000000, 28.5000000000, 28.9375000000, 29.0000000000, 29.0000000000, 29.5000000000, 30.0000000000, 32.0000000000, 33.0000000000, 33.0000000000, 34.0000000000, 34.0000000000, 34.5000000000, 35.0000000000, 35.5000000000, 36.0000000000, 36.0000000000, 36.0000000000, 36.5312500000, 37.0000000000, 37.5001220703, 38.5000000000, 39.0000000000, 40.8750000000, 41.0000000000, 41.0000000000, 42.0000000000, 42.0000305176, 43.0000000000, 43.0000000000, 43.0000000000, 43.9986572266, 44.0000000000, 44.0000000000, 44.0000000000, 44.0000000000, 45.0000000000, 45.0000000000, 46.0000000000, 46.0000000000, 46.0000000000, 46.0002441406, 46.5000000000, 47.0000000000, 48.0000000000, 49.0000000000, 49.0000000000, 50.0000000000, 50.9999389648, 52.0000000000, 53.0000000000, 53.0000000000, 53.9999847412, 54.0000000000, 54.0000152588, 56.9999389648, 57.0000000000, 57.0000000000, 59.0000000000, 60.0000000000, 61.0000000000, 61.0000000000, 61.0000000000, 62.0000000000, 62.9999847412, 63.0000000000, 65.0000000000, 65.0000000000, 65.0000000000, 67.3750000000, 68.0000000000, 68.0000000000, 69.0000305176, 70.0000000000, 71.9999389648, 72.0000000000, 73.0000000000, 74.0000000000, 74.9995117188, 75.0000000000, 76.0000000000, 76.9997558594, 77.0000000000, 77.0000000000, 77.0001144409, 79.0000000000, 79.0000000000, 80.0000000000, 80.0002441406, 85.0000000000, 85.0000000000, 86.0000000000, 87.0000000000, 88.0000000000, 90.9999542236, 91.0000000000, 91.0000000000, 91.0000000000, 92.0000000000, 93.0000305176, 95.0000000000, 96.0000000000, 97.0000000000, 97.0000000000, 97.0000076294, 98.9997558594, 99.0000000000, 100.0000000000, 102.0000000000, 104.0000000000, 105.0000000000, 105.0000000000, 106.0000000000, 106.0000000000, 106.0000000000, 108.9996948242, 109.9989013672, 113.0000305176, 114.0000305176, 116.0000000000, 119.0000000000, 119.0000000000, 121.0000000000, 125.0000000000, 125.9987792969, 126.0000000000, 127.0000000000, 130.0000000000, 135.0000000000, 135.0001220703, 144.9999694824, 146.5000000000, 147.0000000000, 147.0000000000, 150.0000000000, 156.0000000000, 157.0000000000, 165.0000000000, 172.0000000000, 176.9998779297, 178.0000000000, 182.0000000000, 185.0000000000, 187.0000610352, 187.0000610352, 188.9995117188, 189.0000000000, 202.0000000000, 209.0000000000, 227.2500000000, 229.0000000000, 248.0000000000, 271.0000000000, 280.0000000000, 285.0000000000, 298.0000000000, 337.0000000000, 349.0000000000, 351.0000000000]
    double_fault = [0.0000000023, 0.0000000049, 0.0000000049, 0.0000000057, 0.0000000059, 0.0000000061, 0.0000000062, 0.0000000063, 0.0000000063, 0.0000000063, 0.0000000067, 0.0000000069, 0.0000000074, 0.0000000076, 0.0000000077, 0.0000000078, 0.0000000078, 0.0000000083, 0.0000000084, 0.0000000087, 0.0000000089, 0.0000000091, 0.0000000093, 0.0000000094, 0.0000000098, 0.0000000102, 0.0000000107, 0.0000000111, 0.0000000112, 0.0000000113, 0.0000000116, 0.0000000124, 0.0000000130, 0.0000000132, 0.0000000135, 0.0000000138, 0.0000000144, 0.0000000174, 0.0000000181, 0.0000000183, 0.0000000184, 0.0000000195, 0.0000000217, 0.0000000225, 0.0000000227, 0.0000000230, 0.0000000236, 0.0000000240, 0.0000000248, 0.0000000254, 0.0000000256, 0.0000000268, 0.0000000284, 0.0000000284, 0.0000000292, 0.0000000302, 0.0000000302, 0.0000000315, 0.0000000315, 0.0000000317, 0.0000000326, 0.0000000328, 0.0000000329, 0.0000000341, 0.0000000342, 0.0000000372, 0.0000000386, 0.0000000390, 0.0000000402, 0.0000000407, 0.0000000420, 0.0000000423, 0.0000000430, 0.0000000430, 0.0000000432, 0.0000000433, 0.0000000435, 0.0000000454, 0.0000000454, 0.0000000462, 0.0000000482, 0.0000000500, 0.0000000509, 0.0000000528, 0.0000000529, 0.0000000543, 0.0000000547, 0.0000000552, 0.0000000554, 0.0000000566, 0.0000000567, 0.0000000569, 0.0000000593, 0.0000000605, 0.0000000623, 0.0000000641, 0.0000000659, 0.0000000674, 0.0000000683, 0.0000000685, 0.0000000691, 0.0000000692, 0.0000000701, 0.0000000706, 0.0000000724, 0.0000000727, 0.0000000727, 0.0000000732, 0.0000000735, 0.0000000741, 0.0000000748, 0.0000000748, 0.0000000751, 0.0000000752, 0.0000000788, 0.0000000808, 0.0000000821, 0.0000000837, 0.0000000842, 0.0000000846, 0.0000000863, 0.0000000866, 0.0000000867, 0.0000000877, 0.0000000879, 0.0000000884, 0.0000000892, 0.0000000893, 0.0000000900, 0.0000000910, 0.0000000915, 0.0000000936, 0.0000000937, 0.0000000942, 0.0000000974, 0.0000000982, 0.0000000984, 0.0000000987, 0.0000000991, 0.0000000996, 0.0000000998, 0.0000001041, 0.0000001043, 0.0000001052, 0.0000001063, 0.0000001076, 0.0000001121, 0.0000001144, 0.0000001160, 0.0000001169, 0.0000001178, 0.0000001183, 0.0000001191, 0.0000001209, 0.0000001228, 0.0000001242, 0.0000001243, 0.0000001254, 0.0000001264, 0.0000001285, 0.0000001292, 0.0000001295, 0.0000001303, 0.0000001337, 0.0000001355, 0.0000001373, 0.0000001401, 0.0000001414, 0.0000001425, 0.0000001461, 0.0000001488, 0.0000001509, 0.0000001519, 0.0000001522, 0.0000001526, 0.0000001531, 0.0000001554, 0.0000001557, 0.0000001568, 0.0000001597, 0.0000001604, 0.0000001621, 0.0000001631, 0.0000001631, 0.0000001633, 0.0000001636, 0.0000001693, 0.0000001706, 0.0000001711, 0.0000001727, 0.0000001739, 0.0000001754, 0.0000001757, 0.0000001779, 0.0000001839, 0.0000001932, 0.0000001933, 0.0000002007, 0.0000002028, 0.0000002029, 0.0000002034, 0.0000002050, 0.0000002085, 0.0000002096, 0.0000002113, 0.0000002119, 0.0000002125, 0.0000002155, 0.0000002211, 0.0000002260, 0.0000002274, 0.0000002277, 0.0000002280, 0.0000002292, 0.0000002318, 0.0000002363, 0.0000002401, 0.0000002448, 0.0000002464, 0.0000002489, 0.0000002516, 0.0000002537, 0.0000002618, 0.0000002663, 0.0000002684, 0.0000002695, 0.0000002765, 0.0000002775, 0.0000002817, 0.0000002828, 0.0000002873, 0.0000002950, 0.0000002986, 0.0000003020, 0.0000003078, 0.0000003473, 0.0000003568, 0.0000003785, 0.0000003787, 0.0000004206, 0.0000004295, 0.0000004521, 0.0000004548, 0.0000004724, 0.0000004794, 0.0000005010, 0.0000005230, 0.0000005739, 0.0000007003, 0.0000007310]
    long_double_fault = [0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000001, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000002, 0.0000000003, 0.0000000003, 0.0000000003, 0.0000000003, 0.0000000003]
    
class task3:
    double_tlu_time = [0.00000, 200.40000, 401.40000, 1706.00000, 3211.60000, 5415.60000, 8529.40000, 12534.20000, 15941.80000, 23161.40000, 28071.00000, 35895.20000, 47927.20000, 56750.40000, 72893.60000, 91743.80000, 103073.60000, 126235.40000, 168855.00000, 190506.60000, 220085.00000, 258186.00000, 308821.00000, 361761.60000, 416816.20000, 489200.80000, 561392.40000, 632181.00000, 726532.00000, 817874.80000, 947624.40000, 1048529.60000, 1287369.40000, 1334649.20000, 1481941.00000, 1630736.80000, 1787153.20000, 1953030.20000, 2162466.20000, 2295608.00000,]
    double_ldl_time = [0.00000, 200.40000, 605.00000, 901.80000, 1600.80000, 2812.00000, 4711.40000, 6818.00000, 11323.20000, 12733.60000, 15140.00000, 19451.60000, 24064.00000, 29879.20000, 37299.00000, 46222.60000, 52338.80000, 62365.60000, 75199.40000, 90339.80000, 105881.40000, 118416.80000, 133454.40000, 150700.40000, 180197.40000, 204351.20000, 234423.40000, 270920.20000, 305011.00000, 345017.00000, 382517.00000, 425943.40000, 540899.20000, 542048.00000, 613832.00000, 674995.00000, 741682.60000, 817997.00000, 904202.60000, 964063.40000]
    
class task5:
    omega_0_0_iterations = [16, 13, 12, 12, 11, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
    omega_0_5_iterations = [12, 10, 10, 10, 10, 10, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
    omega_1_0_iterations = [7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
    omega_1_5_iterations = [13, 11, 11, 10, 10, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
    omega_2_0_iterations = [16, 13, 12, 11, 11, 11, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
    
    omega_0_8_iterations = [11, 10, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
    omega_1_2_iterations = [11, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
    