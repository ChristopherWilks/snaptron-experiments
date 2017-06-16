#!/usr/bin/env python2.7
import SciServer
from SciServer import Authentication, LoginPortal, Config, CasJobs, SkyQuery, SciDrive, SkyServer
import os
import sys
import json

Authentication_loginName = 'snaptron2_agent';
Authentication_loginPassword = 'snaptron2'
token1 = Authentication.login(Authentication_loginName, Authentication_loginPassword)

CasJobs_TestDatabase = "Snaptron2"
CasJobs_TestQuery = "SELECT top 10 * from SnapCounts"
CasJobs_TestCSVFile = "SciScriptTestFile.csv"

#tables = CasJobs.getTables(context="Snaptron2")
#print tables
df = CasJobs.executeQuery(sql=CasJobs_TestQuery, context=CasJobs_TestDatabase, format="csv")
print df
