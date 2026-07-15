
# --- AESTHETICS ---
SDG_PALETTE <- c(
    "N Africa & W Asia" = "#FB8C00", "S-S Africa"        = "#FFB74D", 
    "N America"         = "#1976D2", "S America"         = "#00BCD4", 
    "C&S Asia"          = "#8BC34A", "E Asia"            = "#55AD31", "SE Asia"           = "#B5C23D",
    "N Europe"          = "#880E4F", "W Europe"          = "#C2185B", "E Europe"          = "#EC407A", "S Europe"          = "#F48FB1",  
    "Oceania"           = "#673AB7"
)

CONT_PALETTE <- c(
    "Africa & W Asia" = "#FB8C00", "N America"       = "#1976D2", "S America"       = "#00BCD4", 
    "Asia"            = "#8BC34A", "Europe"          = "#E91E63", "Oceania"         = "#673AB7", "Global"          = "black"
)

CONT_ORDER <- c("Europe", "N America", "Oceania", "S America", "Asia", "Africa & W Asia")

REGION_ORDER <- c(
    "S Europe", "W Europe", "N Europe", "E Europe", "Oceania", 
    "N America", "S America", "C&S Asia", "E Asia", "SE Asia", 
    "N Africa & W Asia", "S-S Africa"
)

# --- META INFO ---
AGE_WEIGHTS <- c(12000, 10000, 9000, 9000, 8000, 8000, 6000, 6000, 6000, 
                 6000, 5000, 4000, 4000, 3000, 2000, 1000, 500, 500) / 100000

REG_URBAN_2017 <- c(
    "E Europe" = 0.695, "N Europe" = 0.807, "S Europe" = 0.705, "W Europe" = 0.807, 
    "N America" = 0.805, "Oceania" = 0.675, "S America" = 0.839, "C&S Asia" = 0.353, 
    "E Asia" = 0.638, "SE Asia" = 0.492, "N Africa & W Asia" = 0.657, "S-S Africa" = 0.430
)

REG_URBAN_2050 <- c(
    "E Europe" = 0.766, "N Europe" = 0.847, "S Europe" = 0.788, "W Europe" = 0.869, 
    "N America" = 0.866, "Oceania" = 0.727, "S America" = 0.903, "C&S Asia" = 0.470, 
    "E Asia" = 0.896, "SE Asia" = 0.679, "N Africa & W Asia" = 0.700, "S-S Africa" = 0.597
)

HDI_2017 <- c(
    "South Africa" = 0.726,
    "Mauritius" = 0.799,
    "Seychelles" = 0.844,
    "Kenya" = 0.599,
    "Uganda" = 0.557,
    "Kuwait" = 0.838,
    "Turkey" = 0.833,
    "Morocco" = 0.674,
    "Algeria" = 0.746,
    "Israel" = 0.908,
    "Philippines" = 0.694,
    "Singapore" = 0.935,
    "Brunei Darussalam" = 0.837,
    "Thailand" = 0.794,
    "Japan" = 0.920,
    "Republic of Korea" = 0.922,
    "China" = 0.765,
    "Iran (Islamic Republic of)" = 0.793,
    "India" = 0.649,
    "Ecuador" = 0.764,
    "Brazil" = 0.770,
    "Peru" = 0.775,
    "Chile" = 0.862,
    "Costa Rica" = 0.809,
    "Trinidad and Tobago" = 0.806,
    "Uruguay" = 0.827,
    "Argentina" = 0.861,
    "Colombia" = 0.774,
    "Canada" = 0.935,
    "USA" = 0.931,
    "Australia" = 0.940,
    "New Zealand" = 0.936,
    "Belarus" = 0.828,
    "Ukraine" = 0.785,
    "Russian Federation" = 0.838,
    "Czech Republic" = 0.904,
    "Poland" = 0.882,
    "Denmark" = 0.951,
    "Latvia" = 0.870,
    "UK" = 0.937,
    "Estonia" = 0.894,
    "Iceland" = 0.962,
    "Lithuania" = 0.885,
    "Sweden" = 0.949,
    "Norway" = 0.965,
    "Finland" = 0.942,
    "Austria" = 0.922,
    "Ireland" = 0.941,
    "Liechtenstein" = 0.929,
    "Germany" = 0.951,
    "Switzerland" = 0.960,
    "France" = 0.908,
    "The Netherlands" = 0.945,
    "Spain" = 0.902,
    "Croatia" = 0.862,
    "Cyprus" = 0.897,
    "Malta" = 0.901,
    "Italy" = 0.897,
    "Portugal" = 0.864,
    "Slovenia" = 0.918
)