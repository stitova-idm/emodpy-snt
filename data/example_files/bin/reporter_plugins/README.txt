NU_ReportEventRecoder + New Infections Channel in Summary Report
1/25/2020
DanB

NU_ReportEventRecoder - 1/14/2020
This change starts with commit 0dc06d8 from Malaria-Ongoing (Nov 20 2018).
It modifies BaseTextReport.cpp/h such that the method WriteData() opens the
file, writes the data, and closes the file.  This was done to ensure that
the data gets written to the file.  Users were seeing random cases where
ReportEventRecoder did not have any data.

New Infections Channel - 1/25/2020
This change modifies the MalariaSummaryReport by adding a new channel for 
new infections.  The channel is the number of new infections per age bin 
per reporting interval.  

!!!!!!!!!!!!!!!
NOTE: To get data in this channel, you must add NewInfectionEvent to 
Event_Trigger_List in the MalariaSummaryReport definition in custom_reports.json.
The entry should look like:

   "Event_Trigger_List": [
      "EveryUpdate",
      "NewInfectionEvent"
   ],
!!!!!!!!!!!!!!!

Below is some example data where there are three reporting intervals and 11 age bins:

"New Infections": [
    [
        198,
        141,
        85,
        55,
        34,
        36,
        18,
        7,
        4,
        3,
        3
    ],
    [
        482,
        308,
        256,
        153,
        96,
        79,
        41,
        12,
        18,
        11,
        14
    ],
    [
        213,
        149,
        120,
        77,
        41,
        37,
        14,
        12,
        8,
        3,
        10
    ]
]

