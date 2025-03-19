import os
import sys
import random
print(sys.executable)

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from pycirclize import Circos
from pycirclize.utils import fetch_genbank_by_accid
from pycirclize.parser import Genbank
import numpy as np
np.random.seed(0)

import ssl
import urllib.request

from pygenomeviz import GenomeViz

ssl._create_default_https_context = ssl._create_unverified_context

def MakeDatabase(sequences): #Makes the database with all of the sequences.
    directory = r"C:\Users\btb20169\OneDrive - University of Strathclyde\MSci\Fasta Sequences"
    # Dictionary to store SeqRecord objects with strain names as keys
    global Wolbachia
    Wolbachia = {}
    # Loop through each strain and read its sequence
    for strain in sequences:
        print(f"\nReading {strain} Fasta")
        # Parse the sequence for the current strain
        for record in SeqIO.parse(f'{directory}\\{strain}.fasta', "fasta"):
            # Stores the SeqRecord in the dictionary
            Wolbachia[strain] = record

# To access a sequence information, use Wolbachia[strain].x, where x can be id (accession number), description or seq (for sequence).
# E.g. Wolbachia["wMel"].seq gives the wMel sequence
# With a sequence, translate by ammending .translate(to_stop={bool}), .complement(), .count(subsequence), .find(subsequence)
# Certain functions word with the sequence; len(sequence), gc_fraction(sequence), 

def makecontigs(sequence,regions):
    contigs = []
    for region in regions:
        contigs.append(sequence[region[0]:region[1]])
    return contigs

def FastaExport(exportee):
    for strain, record in exportee.items():
        output_directory = r"C:\Users\btb20169\OneDrive - University of Strathclyde\MSci\Fasta Sequences\PhageRegions"
        contigs = []
        for Test in PhageRegions[strain]:
            contigs += PhageRegions[strain][Test]
        # Ensure output directory exists
        os.makedirs(output_directory, exist_ok=True)
        # Generate output file path for this strain with "_Exported" appended
        output_file = os.path.join(output_directory, f"{strain}_Predicted_Prophages.fasta")
        
        records = []  # List to hold SeqRecord objects for this strain
        # Split the sequence at '$' to get individual contigs
        for i, contig in enumerate(contigs):
            if contig.strip():  # Skip empty contigs
                # Generate the header based on whether it's a single contig or multiple
                contig_id = record.id
                # Description format: "Accession Number Description, contig X"
                description = f"{record.description}, contig {i+1}/"+str(len(contigs)) if len(contigs) > 1 else f"{record.description}"
                
                # Create SeqRecord for the contig
                contig_record = SeqRecord(Seq(contig), id=contig_id, description=description)
                records.append(contig_record)

        # Write all records to the output FASTA file for this strain
        with open(output_file, "w") as fasta_out:
            SeqIO.write(records, fasta_out, "fasta")
        
        print(f"\nSequences for {strain} have been exported to {output_file}")

        
def Graph():    
    circos = {}
    gv = GenomeViz(track_align_type="center")
    gv.set_scale_bar()
    for strain in PhageRegions.keys():
        print(f"\nRetrieving {strain} information from Genbank")
        #First, Genbank needs to be downloaded and all regions be annotated as tracks.

        # Download strain genbank using id in the fasta file
        gbk_fetch_data = fetch_genbank_by_accid(Wolbachia[strain].id)
        gbk = Genbank(gbk_fetch_data)        

        #Linear Visualization

        print("Creating linear visualization")
        track = gv.add_feature_track(gbk.name, gbk.get_seqid2size(),align_label=False)
        for seqid, features in gbk.get_seqid2features("CDS").items():
            segment = track.get_segment(seqid)
            segment.add_features(features, plotstyle="bigarrow", fc="limegreen", lw=0.5)
        
        

        # Circos Dictionary
        circos[strain] = Circos(sectors={gbk.name:gbk.range_size})
        circos[strain].text(f"{strain}\n\n{Wolbachia[strain].id}",size=15)
        circos[strain].rect(r_lim=(90, 100), fc="lightgrey", ec="none", alpha=1)
        sector = {}
        sector[strain] = circos[strain].sectors[0]
        
        # Plot forward strand CDS
        print("Plotting forward strand")
        f_cds_track, f_cds_feats = {},{}
        f_cds_track[strain] = sector[strain].add_track((95, 100))
        f_cds_feats[strain] = gbk.extract_features("CDS", target_strand=1)
        f_cds_track[strain].genomic_features(f_cds_feats[strain], plotstyle="arrow", fc="salmon", lw=0.5)

        # Plot reverse strand CDS
        print("Plotting backward strand")
        r_cds_track, r_cds_feats = {},{}
        r_cds_track[strain] = sector[strain].add_track((90, 95))
        r_cds_feats[strain] = gbk.extract_features("CDS", target_strand=-1)
        r_cds_track[strain].genomic_features(r_cds_feats[strain], plotstyle="arrow", fc="skyblue", lw=0.5)

        # Plot 'gene' qualifier label if exists
        print("Labelling genes")
        labels, label_pos_list = [], []
        for feat in gbk.extract_features("CDS"):
            start = int(feat.location.start)
            end = int(feat.location.end)
            label_pos = (start + end) / 2
            gene_name = feat.qualifiers.get("gene", [None])[0]
            if gene_name is not None:
                labels.append(gene_name)
                label_pos_list.append(label_pos)
        f_cds_track[strain].xticks(label_pos_list, labels, label_size=4, label_orientation="vertical")

        # Plot xticks (interval = 100 Kb)
        r_cds_track[strain].xticks_by_interval(
            100000, outer=False, label_formatter=lambda v: f"{v/1000000:.1f} Mb"
        )

        #Next, tracks need to be added for each prediction program used
        minwidth = 80
        maxwidth = 82
        width = 1

        for test in PhageRegions[strain].keys():
            if not test == "Common":
                minwidth += -2.5
                maxwidth += -2.5
                print(f"Plotting {test} phages")
                PhageRegions[strain]["Common"] += PhageRegions[strain][test]
                if test == "PHASTEST":
                    colour = "Violet"
                elif test == "PhageBoost":
                    colour = "Indigo"
                elif test == "PhiSpy":
                    colour = "Blue"
                elif test == "VirSorter":
                    colour = "Green"
                elif test == "VIBRANT":
                    colour = "Yellow"
                elif test == "CheckV":
                    colour = "Orange"
                elif test == "geNomad":
                    colour = "Red"
                elif test == "Known":
                    colour = "#009E73"
                elif test == "WO1-BLAST":
                    colour = (0.5, 0.5, 0.5)
                elif test == "WO2-BLAST":
                    colour = (0.5, 0.5, 0.5)
                elif test == "WO3-BLAST":
                    colour = (0.5, 0.5, 0.5)
                elif test == "WO4-BLAST":
                    colour = (0.5, 0.5, 0.5)
                Phagetracks = {}
                Phagetracks[strain] = {}
                Phagetracks[strain][test] = sector[strain].add_track((minwidth,maxwidth))
                for coords in PhageRegions[strain][test]:
                    if coords[0] > coords[1]:
                        track.add_feature(coords[1], coords[0], width, plotstyle="box", fc=colour)
                        Phagetracks[strain][test].rect(
                            start = float(coords[1]),
                            end = float(coords[0]),
                            fc = colour,
                            alpha = 1,
                            )
                    else:
                        track.add_feature(coords[0], coords[1], width, plotstyle="box", fc=colour)
                        Phagetracks[strain][test].rect(
                            start = float(coords[0]),
                            end = float(coords[1]),
                            fc = colour,
                            alpha = 1,
                            )
        #Goes and plots the common overlaps.
        print("Compiling common phage predictions.")
        minwidth = 85 -2.5
        maxwidth = 87 -2.5
        location_start = 0
        Phagetracks = {}
        Phagetracks[strain] = {}
        Phagetracks[strain]["Common"] = sector[strain].add_track((minwidth,maxwidth))
        limits = []
        overlapping = [0,0]
        for entry in PhageRegions[strain]["Common"]:
            limits.append(entry[0])
            limits.append(entry[1])
        limits.sort()
        for boundary in limits:
            location_end = boundary
            overlaps = 0
            for prediction in PhageRegions[strain]["Common"]:
                if location_start >= prediction[0] and location_end <= prediction[1]:
                    overlaps += 1            
            # Chooses a darker colour if there is more overlap
            if overlaps >= 7:
                colour = (0,0,0)
            elif overlaps == 6:
                colour = (0.15,0.15,0.15)
            elif overlaps == 5:
                colour = (0.3,0.3,0.3)
            elif overlaps == 4:
                colour = (0.45,0.45,0.45)
            elif overlaps == 3:
                colour = (0.6,0.6,0.6)
            elif overlaps == 2:
                colour = (0.75,0.75,0.75)
            elif overlaps == 1:
                colour = (0.9,0.9,0.9)
            else:
                colour = (1,1,1)
                PhageRegions[strain]["Overlaps"] += tuple([overlapping[:]])
                overlapping[0] = location_end
            if overlaps > 0:
                overlapping[1] = location_end
                
                
            track.add_feature(location_start, location_end, width, plotstyle="box", fc=colour)
            Phagetracks[strain]["Common"].rect(
                start = float(location_start),
                end = float(location_end),
                fc = colour,
                alpha = 1,
                )
            location_start = location_end

        
        circos[strain].savefig(f"{strain}_circular.png") #Exported as a PNG file (Good figure for using)
        print(f"{strain}_circular.png created")
        

        
        # Plot xticks (interval = 10 Kb)
        r_cds_track[strain].xticks_by_interval(
            10000, outer=False, label_formatter=lambda v: f"{v/1000:.1f} Kb"
            )

        circos[strain].text(f"{strain}\n\n{Wolbachia[strain].id}",size=300)
        circos[strain].savefig(f"{strain}_circular.svg",figsize=(100, 100)) #Exported as a large SVG file (Lossless; allows for zooming in)

        
    print(PhageRegions[strain]["Overlaps"])
    #Run BLAST aligment and filter?
    
        
    gv.savefig(f"Wolbachias.png")
    print(f"Wolbachias.png created")
        

        
        



PhageRegions = {
    "wAnD":{
        "PHASTEST":([284202,293162],[472597,478452],[1201522,1208586]),
        "PhageBoost":([57594,87939],[239665,280031],[489139,511126],[511673,530917],[1171035,1183271],[1204305,1218431]),
        "PhiSpy":([66384,84222],[284202,290617],[301108,308295],[326066,345139],[472597,472597]),
        "VirSorter": ([281833,327479],[453036,508029],[1198273,1218431]),
        "VIBRANT": ([928821,935121],[1201522,1213507]),
        "CheckV": (),
        "geNomad": ([284416,308943],[0,0]),
        "Known": (),
        "WO1-BLAST": ([287527,289236],[289240,290656],[301108,301703],[326084,327418]),
        "WO2-BLAST": ([68914,69880],[258064,258661],[302711,306299]),
        "WO3-BLAST": ([463452,468038],[0,0]),
        "WO4-BLAST": ([172680,173249],[294074,295558],[476740,479582],[510376,511312],[525409,526820],[525426,526250],[534713,535401],[963560,964526],[1205544,1210123]),
        "Common": (),
        "Overlaps": ()
        },
    "wMel":{
        "PHASTEST":([235113,279207],[548107,573485],[626320,642436]),
        "PhageBoost": ([40152,49668],[62503,74900],[245005,254226],[256859,277369],[549882,564092],[796373,816199]),
        "PhiSpy":([235933,268811],[552599,558081],[575460,582158],[626650,647779]),
        "VirSorter": ([239195,305414],[544371,645205]),
        "VIBRANT": ([239195,278232],[539714,596805]),
        "CheckV": ([240428,280460],[546998,587732]),
        "geNomad": ([243504,273255],[534345,593422],[603458,635557]),
        "Known": (),
        "WO1-BLAST": ([250623,252550],[253867,255241],[255371,256065],[575568,577505],[578204,580534]),
        "WO2-BLAST": ([262946,265307],[266876,268795],[626655,632222],[632302,633890]),
        "WO3-BLAST": ([592576,596442],[0,0]),
        "WO4-BLAST": ([240498,241052],[241174,243818],[488448,493008],[549466,549769],[549797,551023],[551419,553259],[554077,555203],[556401,564091],[567300,568209],[569174,570486],[607321,610281],[610488,612216],[611406,612248],[1245502,1246473]),
        "Common": (),
        "Overlaps": ()
        },
    "wAlbB-HN2016":{
        "PHASTEST":([803774,810191],[1195301,1210760]),
        "PhageBoost":([1189522,1204679],[0,0]),
        "PhiSpy":([1199495,1207674],[0,0]),
        "VirSorter": (),
        "VIBRANT": ([1185734,1212450],[0,0]),
        "CheckV": (),
        "geNomad": (),
        "Known": (),
        "WO1-BLAST": (),
        "WO2-BLAST": (),
        "WO3-BLAST": ([1450291,1454867],[0,0]),
        "WO4-BLAST": ([579209,579948],[587670,588240],[1015421,1016879],[1183597,1184563],[1198644,1201100],[1201795,1202972],[1204790,1205644],[1206827,1207668],[1208017,1208843],[1208813,1212438],[1389561,1390290],[1390406,1391253],[1390435,1392182],[1392193,1395334],[1395449,1396358],[1440156,1441000],[1440157,1441064],[1442985,1447978]),
        "Common": (),
        "Overlaps": ()
        },
    "wAlbB-Uju":{
        "PHASTEST":([185715,196203],[465215,484883],[661345,677790]),
        "PhageBoost":(),
        "PhiSpy":([436746,500246],[664723,677790]),
        "VirSorter": (),
        "VIBRANT": ([665431,686902],[0,0]),
        "CheckV": (),
        "geNomad": ([667124,676797],[0,0]),
        "Known": (),
        "WO1-BLAST": (),
        "WO2-BLAST": ([465392,469859],[470154,471206]),
        "WO3-BLAST": ([681963,686539],[0,0]),
        "WO4-BLAST": ([398583,400041],[441768,442507],[472939,473826],[475166,478307],[478318,480065],[479247,480094],[480210,480939],[668808,670490],[672299,673485],[674177,675672],[676852,679650],[1353421,1354387],[1491784,1492354]),
        "Common": (),
        "Overlaps": ()
        },
    "wCauB":{
        "PHASTEST":([153819, 176444], [653367, 676929], [679646, 704375], [711372, 725769], [1131067, 1154728], [1347897, 1360975], [1410940, 1419296], [1449153, 1471436], [1478471, 1494000], [1520344, 1529957], [1550849, 1575855], [1590988, 1606177], [1626394, 1640802]),
        "PhageBoost": ([20939, 39182], [162604, 205394], [233363, 244934], [668040, 679981], [1057397, 1082053], [1347458, 1354244], [1355286, 1367679], [1412358, 1421903], [1450437, 1465864], [1488404, 1495493], [1516288, 1538604], [1566701, 1573259]),
        "PhiSpy":([153828,205376],[711372,725769],[1393125,1414877],[1452366,1459157],[1524031,1529957],[1554050,1573259]),
        "VirSorter": ([6039,64904],[143081,245804],[373009,392614],[661125,771216],[1122500,1209748],[1404381,1703166]),
        "VIBRANT": ([666654,746553],[1182301,1196524],[1350965,1521287]),
        "CheckV": ([145353,209866],[668583,739232],[1116061,1172067],[1455332,1640802]),
        "geNomad": ([143081,177376],[636099,747388],[1120216,1161147],[1590988,1641068]),
        "Known": (),
        "WO1-BLAST": ([153936,155878],[155881,158214],[711480,716423],[1131184,1133120],[1133123,1135471],[1408234,1410129],[1410948,1412291],[1412425,1413288],[1478579,1480515],[1480559,1483457],[1554158,1556089],[1556448,1557086],[1626502,1631445],[1559620,1560477]),
        "WO2-BLAST": ([718559,719910],[720229,725753],[1137745,1139008],[1355312,1360947],[1486790,1488112],[1488431,1493985],[1563942,1565530],[1565577,1566823],[1567712,1570092],[1571394,1573208],[1633581,1634932],[1635251,1640786]),
        "WO3-BLAST": ([1159255,1160028],[1159346,1159925],[1162019,1173987],[1168309,1168896],[1169482,1170070]),
        "WO4-BLAST": ([150562,149156],[164646,165327],[172354,164647],[178254,177382],[545491,546464],[653401,654227],[653430,655150],[655153,657331],[660288,661044],[661046,661905],[670610,672144],[672856,674538],[677094,679611],[680949,684232],[684230,683768],[685534,686053],[686374,693787],[693790,693099],[700144,698096],[737179,738871],[740282,742478],[742579,743459],[743454,744451],[744972,746535],[1127776,1126291],[1148326,1141483],[1150560,1149717],[1183395,1183914],[1393013,1392445],[1424201,1424720],[1447537,1448424],[1451606,1452116],[1455482,1456824],[1457464,1460568],[1460538,1468218],[1465122,1464436],[1511705,1512451],[1512577,1519287],[1518443,1519279],[1521758,1522645],[1523202,1523781],[1526282,1527624],[1528264,1531368],[1531338,1535922],[1535950,1535239],[1540822,1538683],[1584635,1586852],[1587009,1589652],[1589663,1591869],[1591991,1592841],[1592956,1599257],[1600556,1601075],[1601396,1608809],[1608812,1608121],[1615166,1613118]),
        "Common": (),
        "Overlaps": ()
        },
    "wHa":{
        "PHASTEST":([287518, 308274], [405776, 426796], [428076, 449930]),
        "PhageBoost":([99028, 132243], [265907, 292954], [297406, 309676], [433627, 445136], [542787, 555230], [596532, 610548], [832947, 868544]),
        "PhiSpy":([290026,308274],[405776,420936],[439645,445136]),
        "VirSorter": ([106267,114186],[279189,334825],[375224,515183]),
        "VIBRANT": ([106267,114186],[279189,334825],[375224,515183]),
        "CheckV": ([404977,448905],[0,0]),
        "geNomad": ([391428,447693],[249891,333025]),
        "Known": (),
        "WO1-BLAST": ([290134,292061],[293378,294752],[294882,295576],[418912,416608],[420828,418921]),
        "WO2-BLAST": ([302409,304726],[306339,308258],[411361,407725]),
        "WO3-BLAST": ([246871,250737],[0,0]),
        "WO4-BLAST": ([280789,281631],[280818,281937],[428115,430388],[437961,433625],[441335,438084],[443781,442532],[447872,444473],[451816,449035],[551799,551080],[1272098,1273069]),
        "Common": (),
        "Overlaps": ()
        },
    "wPip":{
        "PHASTEST":([247653, 260088], [320006, 333608], [338019, 360657], [444716, 456479], [467621, 483289], [1375815, 1386259], [1406693, 1416457]),
        "PhageBoost": ([19890,26993],[165191,176501],[293401,321928],[354457,369071],[425387,446644],[484573,506529],[572327,580504],[830067,850284],[1400753,1412009]),
        "PhiSpy":([250652,274242],[320006,360657],[444716,483289],[1376656,1413160]),
        "VirSorter": ([242536,491749],[1371265,1454693]),
        "VIBRANT": ([431883,483288],[1111417,1125407],[1349890,1464203]),
        "CheckV": ([295756,356112],[1377127,1433660]),
        "geNomad": ([248691,367714],[430749,499743],[1361479,1444592]),
        "Known": (),
        "WO1-BLAST": ([258038,255702],[259980,258041],[331558,329222],[333500,331561],[346157,347241],[347244,349580],[396712,397136],[399838,398672],[455177,454320],[456146,455311],[1381658,1380936],[1383137,1381792],[1386151,1384220],[1455157,1454479]),
        "WO2-BLAST": ([251949,250683],[253581,251996],[325570,320021],[327170,325617],[351701,353286],[353333,355788],[357722,360617],[424413,425142],[450286,444737],[451906,450333],[1378734,1377112]),
        "WO3-BLAST": (),
        "WO4-BLAST": ([248536,247685],[272104,274212],[276977,275987],[280140,276972],[284582,280146],[375492,375002],[467659,469927],[470033,470791],[471677,470103],[476035,473054],[479107,476005],[483882,479122],[485575,484660],[491521,490782],[1127652,1126686],[1394392,1396489],[1400723,1401442],[1405308,1400754],[1408378,1405278],[1410707,1409653],[1413391,1411552],[1415392,1414499]),
        "Common": (),
        "Overlaps": ()
        },
    "wRi":{
        "PHASTEST":([581150, 596499], [755998, 761818], [764049, 785272], [1038465, 1048287], [1087609, 1102958], [1337329, 1358163]),
        "PhageBoost": ([126844, 139866], [407581, 427378], [569051, 587684], [724654, 770830], [858946, 880337], [910049, 925318], [977439, 984940], [1059136, 1094143]),
        "PhiSpy":([250652,274242],[320006,360657],[444716,483289],[1376656,1413160]),
        "VirSorter": ([739543,864276],[581150,661408],[1089290,1158934],[1329976,1373865]),
        "VIBRANT": ([755496,793590],[0,0]),
        "CheckV": ([725599,784439],[0,0]),
        "geNomad": ([567035,607370],[739135,785272],[1074268,1113829],[1337329,1362044]),
        "Known": (),
        "WO1-BLAST": ([593748,591424],[596391,594454],[777143,774831],[780166,778226],[1100207,1097883],[1102850,1100913],[1358055,1353279]),
        "WO2-BLAST": ([582807,581155],[587665,583723],[589333,587745],[772776,771262],[1089266,1087614],[1094124,1090182],[1095792,1094204],[1342736,1340271]),
        "WO3-BLAST": ([617570,613704],[741780,742412],[741828,748702],[750693,751466],[1124029,1120163]),
        "WO4-BLAST": ([631410,628451],[633345,631617],[633377,632535],[754245,755147],[760175,763401],[763527,766365],[767282,768787],[768788,768096],[825151,824501],[826689,825202],[1137869,1134910],[1139804,1138076],[1139836,1138994],[1419588,1420559]),
        "Common": (),
        "Overlaps": ()
        },
    "wAnM":{
        "PHASTEST":(),
        "PhageBoost": ([109987,119758],[174224,186797],[364158,385995],[388036,394938],[437836,454361]),
        "PhiSpy":([1009710,1026500],[0,0]),
        "VirSorter": (),
        "VIBRANT": (),
        "CheckV": (),
        "geNomad": (),
        "Known": (),
        "WO1-BLAST": (),
        "WO2-BLAST": ([395327,391829],[396021,395338],[396824,396327]),
        "WO3-BLAST": ([542781,547383],[0,0]),
        "WO4-BLAST": ([132557,130404],[134007,132679],[139152,140118],[429849,431395],[559330,562197]),
        "Common": (),
        "Overlaps": ()
        },
    "wBm":{
        "PHASTEST":(),
        "PhageBoost": (),
        "PhiSpy":([430782,443392],[0,0]),
        "VirSorter": (),
        "VIBRANT": (),
        "CheckV": (),
        "geNomad": (),
        "Known": (),
        "WO1-BLAST": (),
        "WO2-BLAST": (),
        "WO3-BLAST": (),
        "WO4-BLAST": ([507726,508326],[508622,509603]),
        "Common": (),
        "Overlaps": ()
        },
    "wNo":{
        "PHASTEST":([122795,143974],[804982,819223],[1102553,1119411]),
        "PhageBoost": ([106488,128838],[991269,1003648],[1116506,1158363],[1188511,1205397],[1229996,1247681]),
        "PhiSpy":([109939,136893],[808118,819223],[1103588,1149021]),
        "VirSorter": ([111612,170031],[792006,829505],[977439,995610],[1092658,1119913]),
        "VIBRANT": ([806936,820518],[0,0]),
        "CheckV": ([1100293,1167066],[0,0]),
        "geNomad": ([111976,150565],[1102553,1120256]),
        "Known": (),
        "WO1-BLAST": ([134846,132532],[136785,134849],[988270,990198],[991283,993425]),
        "WO2-BLAST": ([127538,123114],[128909,127844],[130583,129225],[813022,814616],[814663,818327]),
        "WO3-BLAST": ([232613,237189],[0,0]),
        "WO4-BLAST": ([140031,141487],[243796,248064],[695519,694553],[808118,809213],[1003124,1003427],[1107529,1108213],[1111880,1107530],[1115238,1112878]),
        "Common": (),
        "Overlaps": ()
        }
    }

def Log():
    with open("log.txt", "w") as log_file:
        for strain in PhageRegions.keys():
            PhageRegions[strain]["Overlaps"]
            log_file.write(strain +" Common Phage Regions:\n")
            for region in PhageRegions[strain]["Overlaps"]:
                log_file.write(str(region[0]) +" - " +str(region[1]) +"\n")
            log_file.write("\n")
    print("Created Log File with overlapping regions")

sequences = PhageRegions.keys()
MakeDatabase(sequences)
input("\nPress enter when Ready to graph.")
Graph()
Log()
     
for Strain in PhageRegions.keys():
    for Test in PhageRegions[Strain].keys():
        PhageRegions[Strain][Test] = makecontigs(Wolbachia[Strain].seq,PhageRegions[Strain][Test])

FastaExport(Wolbachia)



    
