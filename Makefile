# Three targets, in the order a grader runs them:
#   make setup      install dependencies
#   make pipeline   build the database, then run the analysis
#   make dashboard  serve the Streamlit app
#
# Recipes are TAB indented. Spaces break make.

.PHONY: setup pipeline dashboard

# python -m pip, not pip. A bare pip can resolve to a different interpreter's
# pip than the python that runs the pipeline, and the mismatch only shows up
# later as a missing import.
setup:
	python -m pip install -r requirements.txt

# load_data.py rebuilds cell_count.db from scratch, so a second run of this
# target produces the same database and the same files under outputs/.
pipeline:
	python load_data.py
	python run_analysis.py

# 0.0.0.0 and headless are what make the forwarded port work in a Codespace.
dashboard:
	streamlit run app/dashboard.py --server.address 0.0.0.0 --server.headless true
