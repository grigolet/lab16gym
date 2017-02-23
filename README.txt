CONFIG

Place the .txt inside a folder called data_raw-<date>/ and put this folder inside the data/ one.
Next, create a file inside the fits/ folder called data_<date>.json and modify its content (use
data_0220.json as reference).

RUN
Open terminal and type:

python gimmy.py path_to_json_config.json y/n

OPTIONS
y prints the plots to file
n prints only the results (faster)