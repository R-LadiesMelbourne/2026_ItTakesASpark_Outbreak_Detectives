#title: "Outbreak Detectives: Can You Stop the Spread?"
#author: "R-Ladies+ Melbourne"

#{r load-patients}
patients_raw <- load_patients()

#{r look-at-data}
head(patients_raw)

#{r column-names}
names(patients_raw)

#{r check-cleaning}
count_cases(patients)

#{r count-cases}
count_cases(patients)

#{r sickest}
show_sickest(patients)






#{r symptom-table}
symptom_table(patients)

#{r symptom-plot}
symptoms <- symptom_table(patients)

ggplot(symptoms, aes(
  x = n,
  y = reorder(symptom, n)
)) +
  geom_col(fill = "purple") +
  labs(
    title = "Symptoms reported by sick students",
    x = "Number of students",
    y = "Symptom"
  )






#{r temperature}
ggplot(patients, aes(
  x = case_status,
  y = temperature_c,
  fill = case_status
)) +
  geom_boxplot() +
  labs(
    title = "Do sick students have higher temperatures?",
    x = "Student group",
    y = "Temperature (°C)"
  )


#{r How quickly did students become sick - onset}
show_onset(patients)





#{r locations}
locations <- load_locations()

#{r join}
investigation <- combine_evidence(
  patients,
  locations
)

#{r location-table}
location_table(investigation)

#{r location-plot}
locations_summary <- location_table(investigation)

ggplot(locations_summary, aes(
  x = attack_rate,
  y = reorder(location, attack_rate)
)) +
  geom_col(fill = "skyblue") +
  labs(
    title = "Where were sick students exposed?",
    x = "Attack rate",
    y = "Location"
  )








#{r disease-clues}
disease_clues <- load_disease_clues()

disease_clues






#{r tree, fig.width=10, fig.height=7}
show_tree(
  sample_colour = "red",
  reference_colour = "grey30",
  title = "Which reference is closest to the outbreak samples?"
)


 





 
#{r contacts}
contacts <- load_contacts()


#```{r network-groups}
network_groups(contacts)


#{r network, fig.width=10, fig.height=7}
network <- show_network(
  contacts,
  patients,
  sick_colour = "red",
  well_colour = "blue",
  title = "Who had contact with whom?"
)

network
