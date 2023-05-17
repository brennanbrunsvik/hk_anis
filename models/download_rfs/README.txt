BRB2023/02/08. This is a modified version of the code from IRIS EARS to download receiver functions from that database. 
It takes receiver functions from EARS, and puts them in a nice format that works with this HK repository. 
To run: 
- Change your station list in stationList.csv
- Activate a Python environment that has Obspy. 
- Execute: python getEarsRFByStation.py
- Execute: ./stupid_loop.bash


getEarsRFByStation.py downloads data. 
There was some glitch in the code that required re-downloading and processing multiple times. 
I isolated this glitch to the code "stupid_loop.bash", which just re-executes that code many times untill succesfull. 
It's not optimal, but it eventually works! 


Below is the original readme from IRIS EARS

____________
2019-09-12

This bundle contains two Python scripts to facilitate receiver function download from the EARS web site.

Requirements:
    Python 3.7.3 or higher

Event-based receiver functions
------------------------------
Users may download event-based receiver functions using the "Download as ZIP file" link provided on each 
earthquake page, right above the big table, or as an alternative, use the EARS download URL containing the 
event date and time, for example:

http://ears.iris.washington.edu/winkle/earthquakes/2000/01/17/21/18/04

The getEarsRFByEvent.py script provided in this bundle facilitates download of receiver functions for multiple events. 
To use this script:

 1. modify the script and provide the data path

 2. on the EARS main page, search for the desired events using the search form on the
    right panel. select "cvs" as the output option and save the output as "eventList.csv".

 3. run the above python script in the same directory as the eventList.csv file.
    The script will print one pair of "mkdir" and "wget" commands per event to download a zip file of
    all receiver functions .      

 4. create a shell script of the "wget" commands and run the shell script



Station-based receiver functions
---------------------------------
Users may download the station-based receiver functions using the "Download zip of receiver functions as SAC" link 
provided on each station page, just above the table of events, or as an alternative, user may modify the EARS download URL 
containing the network and station codes, for example:

http://ears.iris.washington.edu/receiverFunction.zip?netCode=IU&stacode=SNZO&minPercentMatch=80&gaussian=2.5

The getEarsRFByStation.py script provided in this bundle facilitates download of receiver functions for multiple stations.
To use this script:

 1. modify the script and provide the data path

 2. create a file containing one "network,station" pair per line (you can obtain a complete list of EARS stations using
    "Get a summary of the results for all stations as a CSV text file" link on the bottom panel of the EARS page
    (script expects network and station as the first two values and will ignore the rest of the line) and save the file
    as "stationList.csv".

 3. run the getEarsRFByStation.py python script in the same directory as the stationList.csv file.
    The script will print one pair of "mkdir" and "wget" commands per station to download a zip file of all receiver functions .

 4. create a shell script of the "mkdir" and "wget" commands and run

Questions or comments?
----------------------
If you have any comments or questions, please send an email to product@iris.washington.edu

