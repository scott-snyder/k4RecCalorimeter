#
# File: RecFCCeeCalorimeter/tests/options/ECalEndcapTurbineCaloTool_test.py
# Author: scott snyder <snyder@bnl.gov>
# Date: Feb, 2026
# Purpose: Test for ECalEndcapTurbineCaloTool
#

import os
import Configurables as C

compactFile = 'ALLEGRO_o1_v03.xml'
pathToDetector = os.environ.get('K4GEO','') + '/FCCee/ALLEGRO/compact/' + os.path.splitext(compactFile)[0]

geoSvc = C.GeoSvc('GeoSvc',
                  detectors = [os.path.join (pathToDetector, compactFile)])

ecalEndcapTool = C.ECalEndcapTurbineCaloTool \
    ('ecalEndcapGeometryTool',
     readoutName = 'ECalEndcapTurbine')


appmgr = C.ApplicationMgr \
    (TopAlg = [C.k4__recCalo__ECalEndcapTurbineCaloToolTestAlg
               (ECalEndcapTurbineTool = ecalEndcapTool)],
     ExtSvc = [geoSvc])
