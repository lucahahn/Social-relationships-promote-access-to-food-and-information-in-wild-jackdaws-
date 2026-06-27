
Read me file for: 

Social relationships promote access to food and information in wild jackdaws 
Proceedings of the Royal Society B: Biological Sciences
10.1098/rspb. 2026-1462

Overview
		
	All data and code required to reproduce the analyses presented in the paper are provided in this repository.The script dual_events.R contains 
	all code used for data processing, preparation, analysis, and figures. Parts of the script are included to document the processing 
	steps to derive the datasets for analyses from the original raw data and a larger underlying database. The raw data and full database are 
	not included in this repository. However, all processed datasets required to reproduce the analyses and results reported in the paper are provided.
	To run the analyses, open dual_events.R, navigate to the section entitled "Read datasets for analyses", and then proceed to 
	"(3) STATISTICAL ANALYSES".

	We processed and analysed data in in R (version 4.2.2) via R studio (version 2024.9.1.394).
	
Code file

	"dual_events.R". - contains all code

Data files

	"dual_events.csv" - used for Model 2 (coordination), Model 3 (social tolerance)
	"dual_events_dyads.csv" - used for Model 1 (social preference)
	"visit_data_primary4.csv" - used for Model 4 (displacement)
	"visit_data_primary5.csv" - used for Model 5 (next visitor after dyadic event)
	"visit_data_primary_first_juv2.csv" - used for Model 6 and 7 (social learning)
	"feedercoord.csv" - used for the map 

Variables

	"dual_events.csv"
		initiator - RFID tag of initiator
		joiner - RFID tag of joiner
		position - location/ID of feeding station of dyadic event
		initiator_feeder - perch ID of initiator (.1 = feeder perch; .2 = secondary perch)
		joiner_feeder - perch ID joiner (.1 = feeder perch; .2 = secondary perch)
		start - start of dyadic event
		end - end of dyadic event
		duration - duration of dyadic event (s)
		initiator_bout_start - start of initiator visit
		initiator_bout_end - end of initiator visit
		joiner_bout_start - start of joiner visit
		joiner_bout_end - end of joiner visit
		year - year of data collection
		event_id - unique identifier for events (using primary feeder perch)
		secondary_event_id - unique identifier for events (using secondary feeder perch)
		initiator_perch - perch type of initiator (primary or secondary)
		joiner_perch - perch type of joiner (primary or secondary)
		initiator_JID - individual ID of initiator
		joiner_JID - individual ID of joiner
		primary_perch_JID - individual ID of bird at primary perch
		secondary_perch_JID - individual ID of bird at secondary perch
		initiator_duration - visit duration of initiator
		joiner_duration - visit duration of joiner
		initiator_arriv - arrival of initiator 
		joiner_arriv - arrival of joiner
		initiator_dep - departure of initiator 
		joiner_dep - departure of joiner
		dep_diff - difference in departure time between initiator and joiner (s)
		primary_arriv - arrival of bird at primary perch
		secondary_arriv - arrival of bird at secondary perch
		primary_duration - visit duration of bird at primary perch (s)
		initiator_joiner_ID - dyad ID of initiator (first) and joiner (second)
		dyad_ID - dyad ID of both birds in "alphabetical" order (irrespective of arrival order)
		initiator_partnerID23 - most recent documented pair-bonded partner of initiator in 2023
		joiner_partnerID23 - most recent documented pair-bonded partner of joiner in 2023
		initiator_partnerID24 - most recent documented pair-bonded partner of initiator in 2024
		joiner_partnerID24 - most recent documented pair-bonded partner of joiner in 2024
		initiator_partnerID25 - most recent documented pair-bonded partner of initiator in 2025
		joiner_partnerID25 - most recent documented pair-bonded partner of joiner in 2025
		pair23 - whether or not the dyad was a pair in 2023
		pair24 - whether or not the dyad was a pair in 2024
		pair25 - whether or not the dyad was a pair in 2025
		pair - whether or not the dyad was a pair 
		initiator_year - initiator ID and year of data collection
		joiner_year - joiner ID and year of data collection
		initiator_sex - sex of initiator 
		joiner_sex - sex of joiner 
		sex_combination - sex combination of initiator (first) and joiner (second)
		sex_combination_simple - sex combination of both birds irrespective of arrival order
		initiator_age - age of initiator (years)
		joiner_age - age of joiner (years)
		age_diff - age difference of both birds
		age_diff_abs - age difference of initiator and joiner (absolute value)
		initiator_age_class - age class of initiator (adult or juvenile)
		joiner_age_class - age class of joiner (adult or juvenile)
		age_class_combination - age class combination of both birds 
		age_class_combination_simple - age class combination of both birds irrespective of arrival order
		initiator_motherID - mother of initiator
		joiner_motherID - mother of joiner
		initiator_fatherID - father of initiator 
		joiner_fatherID - father of joiner
		offspring_parent - whether or not initiator was offspring and joiner was a parent
		parent_offspring - whether or not initiator was a parent and joiner was offspring
		parent_offspring_kin - whether or not both birds were parent and offspring
		mother_sibling - whether or not both birds have the same mother
		father_sibling - whether or not both birds have the same father
		sibling - whether or not both birds are siblings
		kin - whether or not both birds are kin (including first-order and second-order kin)
		initiator_box - nest box of initiator in that year 
		joiner_box - nest box of joiner in that year 
		initiator_box_binary - whether or not initiator had a nest box in that year
		joiner_box_binary - whether or not joiner had a nest box in that year
		box_dyad_ID - dyad ID of both nest boxes
		box_distance - distance between both nest boxes 
		neighbour - whether or not two birds where neighbours
		site - study site of dyadic event 
		initiator_pref_site - primary study site of initiator 
		joiner_pref_site - primary study site of joiner
		site_comb - study site combination of initiator and joiner
		site_comb_simple - study site combination of initiator and joiner irrespective of arrival order
		site_comb_binary - whether or not both birds have the same primary study site
		relationship - relationship type between both birds
		social_pedigree - dyadic relatedness between both birds based on a social pedigree
		prev_edgef - association strength (simple ratio index) of both birds at feeders during July of the previous year
		pref_edgef_binary - whether or not both birds associated at feeders during July of the previous year
		day - day of the year
		prop_dual_overlap - relative proportion of dyadic event out of the initiator's total visit time
		primary_duration2 - visit duration of the individual at the primary feeding station that includes a limit for censoring
		primary_duration_cens - identifier of censored values where visit duration measured as 0s
		arriv_diff_z - arrival difference between both birds (z-transformed)
		duration2 - dyadic event duration that includes a limit for censoring
		duration_cens - identifier of censored values where duration measured as 0s

	"dual_events_dyads.csv"
		most variables: see above ("dual_events.csv")
		visit_number_dyad - number of dyadic events for each dyad
		obs - row number (to be used as observation-level random effect)

	"visit_data_primary4.csv" 
		event_id - unique identifier for visit
		feeder	- feeding station ID (only primary perches)
		tag - RFID tag of visiting bird
		start - start of visit
		end - end of visit
		event - event as identified by raw RFID datastream
		time_diff_pre - latency of arrival after previous individual had left (s)
		JID - individual ID of visiting bird
		year - year of data collection
		visit duration - visit duration (s)
		perch - perch type (only primary perches)
		site - study site
		position - location or ID of feeding station 
		day - day of the year
		time - time of day
		hour - hour of the day
		rings - ring combination of visiting bird
		sex - sex of visiting bird
		partnerID - partner of visiting bird
		box - nest box of visiting bird
		motherID - mother of visiting bird 
		fatherID - father of visiting bird
		age - age of visiting bird (years)
		position_year - feeding station location and year of data collection
		pre_JID - ID of bird visiting before current visit
		pre_start - start of previous visit
		pre_end - end of previous visit
		this_start - start of current visit
		this_end - end of current visit
		pre_event_id - unique identifier for previous visit
		next_JID - ID of bird visiting after current visit
		next_start - start of next visit
		next_end - end of next visit
		time_diff_next - latency between end of this and start of next visit
		self_follow - whether or not next visit by the same individual as current visit
		parent - whether subsequent visits are by parent and offspring 
		age_24 - age of individual in 2024 (years)
		JID_position - ID of bird currently visiting and feeding station ID
		JID_position_year - ID of bird currently visiting and feeding station ID and year of data collection
		dual_event_key - whether current visit is a dual event 
		secondary_event_id - unique identifier of current visit at secondary perch 
		duration - duration of dyadic event where applicable 
		arriv_diff - latency of arrival between both birds during dyadic event 
		relationship - relationship type of two birds (solo where no dyadic event)
		secondary_perch_id - ID of secondary perch (dyadic events only)	
		secondary_perch_JID - ID of bird at secondary perch (dyadic events only)
		dyad ID - ID of both birds in case of dyadic event and ID of current and next bird for solo events
		age_class_combination_simple
		pre_queuer_JID - ID of bird queueing at secondary perch before this visit
		pre_relationship - relationship type of birds during dyadic event before this visit (solo = no dyadic event before)
		pre_dyad_ID - dyad ID of birds in dyadic event before this visit
		pre_age_class_combination_simple - age class combination of birds in dyadic event before this visit
		pre_duration - duration of dyadic event before this visit
		pre_arriv_diff - latency of arrival between both birds in dyadic event before this visit
		JID_pre_queuer_JID_same - whether or not individual at feeder perch was a queuer during previous visit
		JID_pre_queuer_JID_same_self_follow - whether or not individual at feeder perch was a queuer during previous visit and whether two subsequent visits of same individual 
		relationship2 - simplified relationship category during this visit (solo, bonded, other)
		pre_relationship2 - simplified relationship category before this visit 
		displaced - whether or not next individual displaced individual currently visiting 
		JID_queuer_JID_displacer_same - whether or not the individual queueing before is now displacing the individual at primary perch
		secondary_perch_age24 - age of bird at secondary perch in 2024
		secondary_perch_age - age of bird at secondary perch in year of data collection (solo if no bird)
		n_contexts - number of contexts that an individual visited in (with bonded partner, with other individual, solo)
		all_contexts - whether or not individual visited in all contexts	
		 
	"visit_data_primary_first_juv2.csv" 
		most variables: see above ("visit_data_primary4.csv")
		day - day of the year at which a juvenile first visited the focal feeding station
		dual_juv - whether or not a juvenile followed its parents at least at one focal feeding station 
		dual_juv_parent - following status of juveniles (2 = "followed", 1 = "sometimes indepenent, 0 = "always independent")
		visit_no - visit number by juvenile at focal feeding station 		
		day_ringed - day of the year that the brood was ringed 
		mother_ID_position_year - mother ID, feeding station position, year of data collection
		father_ID_position_year - father ID, feeding station position, year of data collection
		mother_position_first - day of year that mother first visited focal feeding station
		father_position_first - day of year that father first visited focal feeding station
		mother_days_before - number of days between mother's first visit and juvenile's first visit
		father_days_before - number of days between father's first visit and juvenile's first visit 
		motherID_year - mother ID, year
		fatherID_year - father ID, year
		mother_feeder_use - mother's use of focal feeding station: 0 = never detected in that year, 1 = detected in that year but not at this position, 2 = detected this year and also at this position
		father_feeder_use - father's use of focal feeding station: 0 = never detected in that year, 1 = detected in that year but not at this position, 2 = detected this year and also at this position
		parent_feeder_use - adding numbers from two previous columns together 
		pre_age24 - age of bird visiting before in 2024
		pre_age - age of bird visiting before in year of data collection 
		pre_first_day - first day that the previous visitor was detected at focal feeding station
		pre_days_before - number of days that the previous visitor was first detected before juvenile's first visit
		juv_juv - whether or not bird visiting before and bird currently visiting are juveniles 
		pre_dual_event_key - whether or not event before was a dyadic event 
		day_ringed_z - day of the year that the brood was ringed (z-transformed)

	"feedercoord.csv" 
		position - feeding station ID/location
		latitude - latitude of feeding station location
		longitude - longitude of feeding station location
		feederfun - function of feeding station - dyadic or single-perch
		site - study site 
