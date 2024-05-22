function PeakList = RamanID(S,ConditionName,saveaddress)
% Raman Spectrum Identification
% by Li Lin
% Frb. 5th, 2021
% Arguments:
% S - The data input. 1st column: wavenumbers, 2nd column: intensities
% ConditionName - The name of data
% saveaddress - The address in computer to save results following by "\"
% Outputs:
% PeakList.WaveNumbers - The wavenumbers of identified peaks
% PeakList.Bonds - The meanings of the identified peaks


WN = S(:,1); % Wavenumber Axis
I = S(:,2); % Intensity
s = size(S,1); % data size

MVW = max([round(54335*max(I)^(-1.189)),1]); % moving average width
disp(MVW);
Im = movmean(I,MVW); % moving averaged I
BaseLine = movmean(I,250); % moving averaged I as a base
netI = Im-BaseLine;
stdv = std(netI);
for i = 1:s % remove negative value and baseline noise
    netI(i) = max([netI(i)-stdv*1.5,0]); 
end
NL = 0.005; % noise level
netI_mg = netI/(max(I)-min(I)); % normalize
marker = zeros(s,3);
k = 0;
dI = diff(netI_mg); % derivative
for i = 1:s-2 % find peaks
    if (dI(i) > 0) && (dI(i+1) <= 0) && (netI_mg(i) >= NL)
        k = k+1;
        marker(k,1) = WN(i);
        marker(k,2) = i;
        marker(k,3) = netI_mg(i);
    end
end
limdis = 20; % minimum peak distance [cm-1]
% Remove multiple identifications in one large peak:
if k > 1
    marker = sortrows(marker,3,'descend');
    marker = marker(1:k,:);
    for i = 1:k-1
        if marker(i,1) ~= 0
            for j = i+1:k
                if abs(marker(i,1)-marker(j,1)) <= limdis
                    marker(j,:) = [0,0,0];
                end
            end
        end
    end
    marker = sortrows(marker,3,'descend');
    for i = 1:k
        if marker(i,1) == 0
            break;
        end
    end
    k = i-1;
end
marker = marker(1:k,:);
PeakList.Wavenumbers = marker(:,1);
if k == 0
    PeakList.Bonds = 'No peaks found.';
else
    lb = boundarylib();
    meaning = dictionary();
    t = 0;
    bondlist = cell(size(lb,1),1);
    libnum = zeros(size(lb,1),1);
    for i = 1:size(marker,1)
        for j = 1:size(lb,1)
           if (marker(i,1) >= lb(j,1)) && (marker(i,1) <= lb(j,2))
               t = t+1;
               libnum(t) = j;
               bondlist{t} = meaning{j};
               break;
           end
        end
        if t ~= i
            t = t+1;
            bondlist{t} = 'Not in the library';
        end
    end
    bondlist = bondlist(1:t);
    libnum = libnum(1:t);
    PeakList.Bonds = bondlist;
end
close all;
figure('position',[10,200,1200,600]);
plot(WN,I,'k','LineWidth',0.5);
hold on;
vloc = flipud((linspace(min(I)+0.01*max(I),max(I)-0.01*max(I),size(marker,1)))');
for i = 1:size(marker,1)
    xline(marker(i,1),'r','LineWidth',1);
    text(marker(i,1),vloc(i),...
        strcat({'  '},num2str(marker(i,1)),{'   '},bondlist{i}),...
        'color','r','FontSize',10,...
        'FontWeight','bold');
end
hold off;
grid on;
xlabel('Wavenumber [cm^-^1]');
ylabel('Intensity [a.u.]');
xlim([min(WN),max(WN)]);
title(ConditionName);
save(strcat(saveaddress,ConditionName,'-Results.mat'),'PeakList');
saveas(gcf,strcat(saveaddress,ConditionName,'-Results.fig'));
saveas(gcf,strcat(saveaddress,ConditionName,'-Results.jpg'));
disp('Complete');
end

function meaning = dictionary()
meaning = [...
"	Lattice vibrations in crystals, LA modes 	"
"	\nu(XmetalO) 	"
"	\delta(CC) aliphatic chains	"
"	\nu(SeSe) 	"
"	\nu(SS) 	"
"	\nu(SiOSi) 	"
"	\nu(CI) 	"
"	SS bridges	"
"	\nu(CBr) 	"
"	SS stretch 	"
"	\nu(CCl)	"
"	\nu(CC) alicyclic, aliphatic chain vibrations	"
"	Phenylalanine 	"
"	\nu(CS) aliphatic 	"
"	Tyrosine	"
"	Tyrosine 	"
"	CS stretch 	"
"	CS str.	"
"	\nu (\delta (CCN), Vinyl & Porphyrin	"
"	CN str.	"
"	Phospholipid (choline)[26-28]	"
"	Tryptophan 	"
"	Nucleic acids, Trp	"
"	CH2 Rock, Sym. breathing	"
"	\nu(COC)	"
"	tyr	"
"	Tyrosine 	"
"	Phosphodiester BkB 	"
"	Tyrosine	"
"	Tyr, proline, glycogen [16]	"
"	\nu(OO)	"
"	Tyrosine 	"
"	tyr	"
"	Phosphate group	"
"	Most probably the amino acids, polysaccharides collagen	"
"	Tyr, Lipids,Carbohydrates,Collagen [26] C-C-N+, C-O-C ring, C-C	"
"	Phosphodiester BkB, deoxyribose 	"
"	C-C stretch of proline, glucose, lactic acid [16]	"
"	C\alphaC stretch 	"
"	V (CC), ahelix 	"
"	Hydroxyproline,Collagen backbone [27-28] CH=CH bending	"
"	?C-H out of plane deformation C-C Asym. Str.	"
"	C-C symmetric stretching, glucose-I-phosphate, sym. breathing mode of phenylalanine [16]	"
"	Phenylalanine 	"
"	\nu(C=S) 	"
"	\nu(CC) aromatic ring chain vibrations	"
"	Phenylalanine	"
"	Phe	"
"	Symmetric CC aromatic ring breathing	"
"	Phenylalanine [13, 29, 30] sym. ring breathing of protein [31]	"
"	Tryptophan 	"
"	Phenylalanine 	"
"	Phe	"
"	Collagen	"
"	CC str.	"
"	\nu(COC) asym	"
"	Lipids,Collagen [26, 27] C-C str.	"
"	V (CO),  V (CC) 	"
"	Triglycerides	"
"	\nu(CS) aromatic 	"
"	Proteins (C-C,C-N str.)[26, 32, 33], P=O sym. from nucleic acids and phospholipids	"
"	CC stretch, CC skeletal stretch trans, PO2 symmetric	"
"	Symmetric PO2- stretching vibration of the DNA	"
"	Phospholipids, O-P-O sym. str.[27], P=O sym. from nucleic acids,cell membrane phospholipids	"
"	PO2	"
"	V (CN) BkB V (CC) Lipid side chain 	"
"	CN, CC str.	"
"	C-C stretching, trans	"
"	Palmitic acid, fatty acid	"
"	RibosePhosphate 	"
"	CC (&CN) stretching of proteins glycogen, carotenoids most likely a cellular pigment	"
"	C?C stretch	"
"	RibosePhosphate 	"
"	CC, CN str.	"
"	L-Tryptophan [31]	"
"	Tyrosine, Phenylalanine	"
"	C-H in-plane bending	"
"	Tyrosine	"
"	Cytosine, Guanine	"
"	Thymine (T), Guanine (G), Cytosine (C), Phenylalanine 	"
"	Antisymmetric phosphate vibrations	"
"	C-C6H5 Phe, Trp [26]	"
"	Phenylalanine, Tyrosine 	"
"	Amide III	"
"	\betasheet	"
"	Cytosine (C) 	"
"	Amide III	"
"	Phospholipid, O-P-O antisym. Stretch [8] Amide III \betasheet [16]	"
"	Antisymmetric phosphate PO2- (antisymmetric) stretching modes (nucleic acids typical for malignant tissues), the PO2? groups of phospholipids do not contribute to these bands, Amide III (b-sheet and random coils)	"
"	Amide III \beta	"
"	Random coil 	"
"	Nucleic acids (Try, Ala),Proteins (Amide III \beta sheet or random coil), Lipid, phospholipid =C-H bend [26, 27] 	"
"	Amide III \nu (CN) and \delta (NH) of proteins	"
"	Amide III, \delta(C=CH2) 	"
"	Fatty acids, =C-H bend [26, 27]	"
"	?helix 	"
"	Amide III [16] \alpha helix, P=O asymmetric stretch due to nucleic acids	"
"	Amide III \alpha 	"
"	CH twist and bend	"
"	 \delta(CH2) 	"
"	Amide III band 40% CN stretch, 30% NH bend	"
"	CH2 deformation of lipids	"
"	Amide III, \delta (N-H)-30%, \alpha-helix, \nu (C-N)-40% & \delta(CH3)	"
"	Lipids, phospholipids [27]C-H2 twist, collagen, protein amide III [16], DNA [16]	"
"	CH3,CH2 twisting, wagging &,or bending mode of collagens & lipids	"
"	Guanine (G) 	"
"	CH2 Deformation	"
"	Trp, Ca-H def	"
"	Tryptophan 	"
"	\nu(C(NO2)) 	"
"	Tryptophan C?H 	"
"	Adenine (A) 	"
"	CH residual vibrations	"
"	CH3-(C?O),	"
"	Trp, Ca-H def	"
"	Thymine (T), Adenine (A), Guanine (G) 	"
"	CH3 in-phase deformation	"
"	\delta(CH3) 	"
"	CH str.	"
"	\delta(CH2) \delta(CH3) asym 	"
"	\nu(N=N) aromatic 	"
"	\delta (CN) bending, \delta(CH)3 out-of-phase deformation	"
"	Fatty acids, triglycerides, CH2 or CH3 deformations [33]	"
"	CH2 bending mode in normal tissue	"
"	CH def 	"
"	\delta(CH) proteins, \delta(CH) Lipids 	"
"	\nu(CC) aromatic ring chain vibrations	"
"	CH2 bending mode in malignant tissues, bending modes of methyl groups (vibrational modes of collagen)	"
"	Proteins [30, 31] C-H wag, CH2 or CH3 def., Phospholipids, CH2 scissoring [34]	"
"	CH def.,bending	"
"	\nu(CC) aromatic ring chain vibrations	"
"	\nu(C=C) 	"
"	Betacarotene CC stretching mode	"
"	Amide II, Shift to 1548, (C?C) stretch	"
"	\nu(C(NO2)) asym 	"
"	NH and NH2 in cytosine, cytidine	"
"	Amide II, in plane \delta (N-H) bending: 60%; \nu (C-N):40%;	"
"	Amide II band 60% NH bend and 40% CN stretch	"
"	Parallel , Antiparallel \beta sheet structure 	"
"	\nu(N=N) aliphatic 	"
"	Amide II, proteins [27, 33], amide II ?sheet [35]	"
"	Hemoglobin 	"
"	Amide II	"
"	Adenine (A), Guanine (G) 	"
"	Adenine (A), Guanine (G), Purine 	"
"	\nu(CC) aromatic ring chain vibrations	"
"	Tryptophan	"
"	Phenylalanine, Tryptophan 	"
"	Phosphorylated amino acids and proteins	"
"	Phenylalanine, Tryptophan 	"
"	Amide II [16], aromatic amino acids within proteins,[33] nucleic acids [26, 27]	"
"	C-C stretching, C?H bending	"
"	Tyrosine	"
"	\nu(CC) aromatic ring chain vibrations	"
"	Adenine (A), Cytosine (G) 	"
"	Phenylalanine	"
"	Phenylalanine 	"
"	CO stretching, C?C bending	"
"	Tryptophan, Phenylalanine, Tyrosine 	"
"	\nu(C=N) 	"
"	Antiparallel \beta sheet 	"
"	Tyrosine, Tryptophan, C=C (protein)	"
"	Tryptophan	"
"	Parallel \beta sheet	"
"	Amide I in \alpha-helix	"
"	Unordered	"
"	\delta(H2O) 	"
"	Amide I 	"
"	Amide I band 80% C=O stretch	"
"	?helix 	"
"	Amide I \nu (C=O) of proteins	"
"	Turn 	"
"	Amide I 	"
"	Unsaturated fatty acids, triglycerides (C=C) str.[33], Amide I \alpha helix	"
"	Amide I vibration mode of structural proteins C=C cis, lipids, fatty acids	"
"	Amide I, \beta-sheet, \nu (C?O) 80%	"
"	Proteins, Amide I ?sheet, cholesterol esters [16]	"
"	Proteins Amide I turn [23] , Unsaturated fatty acids [26, 31], (C=O) str., (C-H) def.,(C=C) str.[26, 31], collagen, elastin [27]	"
"	Weak Antiparallel \beta sheet 	"
"	Turn 	"
"	\nu(C=O) 	"
"	\nu (C?C)	"
"	(C=O) stretching, triglycerides [33]	"
"	\nu(C\congC) 	"
"	\nu(C\congN) 	"
"	\nu(SH) 	"
"	1378??cm-1 bend overtone	"
"	\nu(C-H) 	"
"	Fatty acids, triglycerides, C-H2 sym. str.	"
"	\nu(CH2)	"
"	Fatty acids, triglycerides, C-H2 sym. str.	"
"	Lipids [31], C-H2 antisym. str.	"
"	\nu (CH2, FR)	"
"	CH stretch of lipids and proteins	"
"	CH band of lipids and proteins	"
"	Proteins,Lipids, CH3 sym. str.[27, 31]	"
"	\nu (CH3, FR)	"
"	Proteins,Lipids, CH3 sym. str.[27, 31]	"
"	\nu(=(CH)) 	"
"	=C-H, lipids, fatty acids	"
"	Lipids [27, 31] =C-H str.	"
"	CH3-(C?O)	"
"	Nucleic acids,Proteins [31] C-H aromatic	"
"	\nu(OH) 	"
"	\nu(O-H) water band	"
"	Histidine	"
"	\nu(O?H) water band	"
"	O-H, Liquid water	"
"	\nu(\cong(CH)) 	"
"	\nu(NH) 	"
"	O-H, Liquid water	"];
end

function lb = boundarylib()
lb = [...
10.00	200.00
150.00	450.00
250.00	400.00
290.00	330.00
430.00	550.00
450.00	550.00
480.00	660.00
500.00	550.00
500.00	700.00
508.00	508.00
550.00	800.00
600.00	1300.00
620.00	620.00
630.00	790.00
640.00	640.00
643.00	643.00
667.00	667.00
672.00	672.00
676.00	676.00
712.00	712.00
721.00	721.00
750.00	750.00
751.00	751.00
754.00	754.00
800.00	970.00
820.00	820.00
830.00	830.00
830.00	830.00
831.00	831.00
840.00	840.00
845.00	900.00
850.00	850.00
852.00	852.00
860.00	860.00
870.00	870.00
883.00	883.00
895.00	895.00
917.00	917.00
939.00	939.00
940.00	940.00
958.00	958.00
973.00	973.00
997.00	997.00
1000.00	1000.00
1000.00	1250.00
1000.00	1000.00
1004.00	1004.00
1004.00	1004.00
1004.00	1004.00
1004.00	1004.00
1011.00	1011.00
1030.00	1030.00
1030.00	1030.00
1039.00	1039.00
1053.00	1053.00
1060.00	1150.00
1064.00	1064.00
1065.00	1065.00
1073.00	1073.00
1080.00	1100.00
1080.00	1158.00
1088.00	1088.00
1091.00	1091.00
1091.00	1091.00
1097.00	1097.00
1126.00	1126.00
1126.00	1126.00
1128.00	1128.00
1130.00	1130.00
1144.00	1144.00
1155.00	1155.00
1156.00	1156.00
1157.00	1157.00
1157.00	1157.00
1160.00	1160.00
1170.00	1200.00
1173.00	1173.00
1174.00	1174.00
1175.00	1175.00
1176.00	1176.00
1187.00	1187.00
1189.00	1189.00
1205.00	1205.00
1214.00	1214.00
1229.00	1235.00
1230.00	1230.00
1230.00	1230.00
1238.00	1238.00
1241.00	1241.00
1242.00	1242.00
1243.00	1253.00
1248.00	1248.00
1260.00	1260.00
1266.00	1266.00
1267.00	1267.00
1270.00	1300.00
1276.00	1276.00
1278.00	1278.00
1280.00	1280.00
1295.00	1295.00
1300.00	1300.00
1301.00	1301.00
1301.00	1301.00
1304.00	1304.00
1310.00	1310.00
1320.00	1320.00
1338.00	1338.00
1339.00	1339.00
1340.00	1360.00
1340.00	1380.00
1341.00	1341.00
1342.00	1342.00
1347.00	1347.00
1358.00	1358.00
1370.00	1370.00
1376.00	1376.00
1378.00	1378.00
1380.00	1380.00
1388.00	1388.00
1400.00	1470.00
1410.00	1440.00
1428.00	1471.00
1437.00	1444.00
1438.00	1438.00
1447.00	1447.00
1448.00	1448.00
1450.00	1450.00
1452.00	1452.00
1453.00	1453.00
1471.00	1471.00
1500.00	1500.00
1500.00	1900.00
1516.00	1516.00
1527.00	1527.00
1530.00	1590.00
1540.00	1540.00
1548.00	1548.00
1550.00	1550.00
1550.00	1550.00
1550.00	1580.00
1558.00	1558.00
1566.00	1566.00
1566.00	1566.00
1578.00	1578.00
1578.00	1578.00
1580.00	1580.00
1582.00	1582.00
1582.00	1582.00
1583.00	1583.00
1584.00	1584.00
1584.00	1584.00
1587.00	1587.00
1590.00	1590.00
1600.00	1600.00
1603.00	1603.00
1605.00	1605.00
1605.00	1605.00
1605.00	1605.00
1610.00	1616.00
1610.00	1680.00
1612.00	1640.00
1615.00	1615.00
1618.00	1621.00
1626.00	1640.00
1639.00	1639.00
1640.00	1651.00
1640.00	1640.00
1643.00	1643.00
1650.00	1650.00
1650.00	1657.00
1652.00	1652.00
1655.00	1675.00
1656.00	1656.00
1658.00	1658.00
1661.00	1661.00
1667.00	1667.00
1667.00	1680.00
1667.00	1680.00
1670.00	1690.00
1680.00	1696.00
1680.00	1820.00
1732.00	1732.00
1732.00	1732.00
2100.00	2250.00
2220.00	2255.00
2550.00	2600.00
2727.00	2727.00
2800.00	3000.00
2845.00	2845.00
2850.00	2850.00
2854.00	2854.00
2888.00	2888.00
2891.00	2891.00
2905.00	2905.00
2911.00	2911.00
2931.00	2931.00
2934.00	2934.00
2940.00	2940.00
3000.00	3100.00
3009.00	3009.00
3009.00	3009.00
3060.00	3060.00
3067.00	3067.00
3100.00	3650.00
3104.00	3104.00
3110.00	3160.00
3156.00	3156.00
3288.00	3288.00
3300.00	3300.00
3300.00	3500.00
3444.00	3444.00
];
end