#
# File: RecCalorimeter/tests/options/CaloCellConstantsSvc_test.py
# Author: scott snyder <snyder@bnl.gov>
# Date: Jan, 2026
# Purpose: Test for CaloCellConstantsSvc
#

import os, sys
for p in os.environ['LD_LIBRARY_PATH'].split(':'):
    print('--- path ' + p, flush=True)
    os.system ('ls ' + p)
    pp = os.path.join (p, 'k4RecCalorimeter.confdb')
    if os.path.exists (pp):
        print('--- file ' + pp, flush=True)
        os.system ('cat ' + pp)
    pp = os.path.join (p, 'k4RecCalorimeter.components')
    if os.path.exists (pp):
        print('--- file ' + pp, flush=True)
        os.system ('cat ' + pp)

import Configurables as C
appmgr = C.ApplicationMgr()
appmgr.TopAlg += [C.k4__recCalo__CaloCellConstantsSvcTestAlg()]
sys.exit(1)
