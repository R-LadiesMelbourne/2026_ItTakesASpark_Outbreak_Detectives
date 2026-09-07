
# ============================================================
# Outbreak Detectives - helper functions
# R-Ladies+ Melbourne
# ============================================================

load_patients <- function() {
  readr::read_csv(
    "case_files/mystery_A/01_patient_records_messy.csv",
    show_col_types = FALSE
  )
}

load_locations <- function() {
  readr::read_csv(
    "case_files/mystery_A/02_locations.csv",
    show_col_types = FALSE
  )
}

load_contacts <- function() {
  readr::read_csv(
    "case_files/mystery_A/04_contacts.csv",
    show_col_types = FALSE
  )
}

load_disease_clues <- function() {
  readr::read_csv(
    "disease_reference.csv",
    show_col_types = FALSE
  )
}

clean_data <- function(data) {

  yes_no_columns <- c(
    "hospitalised", "fever", "cough", "sore_throat", "headache",
    "muscle_aches", "vomiting", "diarrhoea", "loss_of_smell",
    "bleeding", "confusion"
  )

  # Include strawberry_tongue if it exists in the supplied data.
  if ("strawberry_tongue" %in% names(data)) {
    yes_no_columns <- c(yes_no_columns, "strawberry_tongue")
  }

  data |>
    dplyr::mutate(
      sex = dplyr::case_when(
        stringr::str_to_lower(sex) %in% c("f", "female") ~ "Female",
        stringr::str_to_lower(sex) %in% c("m", "male") ~ "Male",
        TRUE ~ "Another identity / not recorded"
      ),
      case_status = dplyr::if_else(
        stringr::str_to_lower(case_status) %in% c("sick", "yes"),
        "Sick",
        "Well"
      ),
      dplyr::across(
        dplyr::all_of(yes_no_columns),
        ~ dplyr::if_else(
          stringr::str_to_lower(.x) %in% c("yes", "y"),
          "Yes",
          "No"
        )
      )
    )
}

count_cases <- function(data) {
  data |>
    dplyr::count(case_status, name = "students")
}

show_sickest <- function(data, number = 6) {
  data |>
    dplyr::filter(case_status == "Sick") |>
    dplyr::arrange(dplyr::desc(temperature_c)) |>
    dplyr::select(patient_id, first_name, temperature_c, hospitalised) |>
    utils::head(number)
}

symptom_table <- function(data) {

  symptoms <- c(
    "fever", "cough", "sore_throat", "headache",
    "muscle_aches", "vomiting", "diarrhoea",
    "loss_of_smell", "bleeding", "confusion"
  )

  if ("strawberry_tongue" %in% names(data)) {
    symptoms <- c(symptoms, "strawberry_tongue")
  }

  data |>
    dplyr::filter(case_status == "Sick") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(symptoms),
      names_to = "symptom",
      values_to = "present"
    ) |>
    dplyr::filter(present == "Yes") |>
    dplyr::count(symptom, sort = TRUE)
}

show_symptoms <- function(data) {

  symptom_table(data) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = n,
        y = forcats::fct_reorder(symptom, n)
      )
    ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = "What symptoms are most common?",
      x = "Number of sick students",
      y = NULL
    ) +
    ggplot2::theme_minimal()
}

show_temperature <- function(data) {

  data |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = case_status,
        y = temperature_c,
        fill = case_status
      )
    ) +
    ggplot2::geom_boxplot(show.legend = FALSE) +
    ggplot2::labs(
      title = "Do sick students have higher temperatures?",
      x = NULL,
      y = "Temperature (°C)"
    ) +
    ggplot2::theme_minimal()
}

show_onset <- function(data) {

  data |>
    dplyr::filter(case_status == "Sick") |>
    dplyr::summarise(
      sick_students = dplyr::n(),
      typical_temperature = round(stats::median(temperature_c, na.rm = TRUE), 1),
      typical_hours_to_symptoms = round(stats::median(hours_to_symptoms, na.rm = TRUE), 1)
    )
}


combine_evidence <- function(patients, locations) {
  dplyr::left_join(patients, locations, by = "patient_id")
}

location_table <- function(data) {

  data |>
    dplyr::group_by(location) |>
    dplyr::summarise(
      students = dplyr::n(),
      sick = sum(case_status == "Sick"),
      attack_rate = sick / students,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(attack_rate), dplyr::desc(sick))
}

show_locations <- function(data) {

  location_table(data) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = attack_rate,
        y = forcats::fct_reorder(location, attack_rate)
      )
    ) +
    ggplot2::geom_col() +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::labs(
      title = "Where were sick students exposed?",
      x = "Attack rate",
      y = NULL
    ) +
    ggplot2::theme_minimal()
}


show_tree <- function(
  sample_colour = "red",
  reference_colour = "grey30",
  title = "Which reference is closest to the outbreak samples?"
) {

  alignment_file <- "case_files/mystery_A/03_sequencing_alignment.fasta"

  dna <- ape::read.dna(
    alignment_file,
    format = "fasta"
  )

  genetic_distance <- ape::dist.dna(
    dna,
    model = "TN93",
    pairwise.deletion = TRUE
  )

  tree <- ape::nj(genetic_distance)
  tree <- phangorn::midpoint(tree)

  set.seed(2026)

  support <- ape::boot.phylo(
    phy = tree,
    x = dna,
    FUN = function(x) {
      phangorn::midpoint(
        ape::nj(
          ape::dist.dna(x, model = "TN93")
        )
      )
    },
    B = 200,
    quiet = TRUE
  )

  tip_colours <- ifelse(
    grepl("^REF_", tree$tip.label),
    reference_colour,
    sample_colour
  )

  graphics::plot(
    tree,
    type = "phylogram",
    direction = "rightwards",
    tip.color = tip_colours,
    cex = 0.8,
    label.offset = 0.002,
    main = title
  )

  ape::nodelabels(
    text = ifelse(support >= 70, support, ""),
    frame = "none",
    cex = 0.65
  )

  graphics::legend(
    "topleft",
    legend = c("Student sample", "Reference sample"),
    col = c(sample_colour, reference_colour),
    pch = 19,
    bty = "n"
  )

  invisible(tree)
}


show_network <- function(
  contacts,
  patients,
  sick_colour = "red",
  well_colour = "blue",
  title = "Who had contact with whom?"
) {

  graph <- igraph::graph_from_data_frame(
    contacts,
    directed = FALSE,
    vertices = patients |>
      dplyr::select(patient_id, first_name, case_status)
  )

  comp <- igraph::components(graph)$membership

  main_nodes <- which(comp == comp["P01"])
  food_nodes <- which(comp == comp["P29"])
  isolated_nodes <- which(igraph::V(graph)$name == "P36")

  set.seed(2026)

  main_layout <- igraph::norm_coords(
    igraph::layout_with_fr(
      igraph::induced_subgraph(graph, main_nodes)
    ),
    xmin = -1, xmax = 1,
    ymin = -1, ymax = 1
  )

  layout <- matrix(
    NA,
    nrow = igraph::vcount(graph),
    ncol = 2
  )

  layout[main_nodes, ] <- main_layout

  layout[food_nodes, ] <- matrix(
    c(
      1.65, -0.65,
      1.65, -1.00
    ),
    ncol = 2,
    byrow = TRUE
  )

  layout[isolated_nodes, ] <- c(-1.55, 0.8)

  igraph::V(graph)$display_label <- igraph::V(graph)$name

  # P01 is the first recognised school case in the Strep A version.
  igraph::V(graph)$display_label[
    igraph::V(graph)$name == "P01"
  ] <- "P01\nFirst reported case"

  igraph::V(graph)$node_size <- ifelse(
    igraph::V(graph)$name == "P01",
    8,
    5
  )

  ggraph::ggraph(
    graph,
    layout = "manual",
    x = layout[, 1],
    y = layout[, 2]
  ) +
    ggraph::geom_edge_link(
      colour = "grey75",
      linewidth = 0.7
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(
        colour = case_status,
        size = node_size
      )
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = display_label),
      repel = TRUE,
      size = 3.5
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Sick" = sick_colour,
        "Well" = well_colour
      ),
      name = NULL,
      labels = c(
        "Sick student",
        "Well student"
      )
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::labs(
      title = title
    ) +
    ggraph::theme_graph() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "left"
    )
}

network_groups <- function(contacts) {
  graph <- igraph::graph_from_data_frame(
    contacts,
    directed = FALSE
  )

  sort(
    igraph::components(graph)$csize,
    decreasing = TRUE
  )
}
