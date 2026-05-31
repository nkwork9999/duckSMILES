pub(crate) struct NamedSmiles {
    pub(crate) name: &'static str,
    pub(crate) smiles: &'static str,
}

pub(crate) const METHANE: &str = "C";
pub(crate) const WATER: &str = "O";
pub(crate) const ETHANOL: &str = "CCO";
pub(crate) const BENZENE: &str = "c1ccccc1";
pub(crate) const PHENOL: &str = "c1ccccc1O";
pub(crate) const ASPIRIN: &str = "CC(=O)Oc1ccccc1C(=O)O";
pub(crate) const SALICYLIC_ACID: &str = "OC(=O)c1ccccc1O";
pub(crate) const PARACETAMOL: &str = "CC(=O)Nc1ccc(O)cc1";
pub(crate) const IBUPROFEN: &str = "CC(C)Cc1ccc(C(C)C(=O)O)cc1";
pub(crate) const CAFFEINE: &str = "Cn1c(=O)c2c(ncn2C)n(C)c1=O";
pub(crate) const CHOLESTEROL: &str = "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C";

pub(crate) const ASPIRIN_SIMILARITY_SET: &[NamedSmiles] = &[
    NamedSmiles {
        name: "aspirin",
        smiles: ASPIRIN,
    },
    NamedSmiles {
        name: "salicylic acid",
        smiles: SALICYLIC_ACID,
    },
    NamedSmiles {
        name: "paracetamol",
        smiles: PARACETAMOL,
    },
    NamedSmiles {
        name: "ibuprofen",
        smiles: IBUPROFEN,
    },
    NamedSmiles {
        name: "caffeine",
        smiles: CAFFEINE,
    },
    NamedSmiles {
        name: "benzene",
        smiles: BENZENE,
    },
    NamedSmiles {
        name: "methane",
        smiles: METHANE,
    },
];

pub(crate) const ECFP_REPORT_SET: &[NamedSmiles] = &[
    NamedSmiles {
        name: "methane",
        smiles: METHANE,
    },
    NamedSmiles {
        name: "water",
        smiles: WATER,
    },
    NamedSmiles {
        name: "ethanol",
        smiles: ETHANOL,
    },
    NamedSmiles {
        name: "benzene",
        smiles: BENZENE,
    },
    NamedSmiles {
        name: "phenol",
        smiles: PHENOL,
    },
    NamedSmiles {
        name: "aspirin",
        smiles: ASPIRIN,
    },
    NamedSmiles {
        name: "caffeine",
        smiles: CAFFEINE,
    },
    NamedSmiles {
        name: "cholesterol",
        smiles: CHOLESTEROL,
    },
];
