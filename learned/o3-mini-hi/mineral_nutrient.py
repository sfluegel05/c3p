"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: Mineral nutrient
Definition: A mineral that is an inorganic nutrient which must be ingested and absorbed
in adequate amounts to satisfy a wide range of essential metabolic and/or structural functions
in the human body.
Examples include various salts containing metal ions such as calcium, iron, magnesium, etc.
"""

from rdkit import Chem

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    The heuristic here is that many mineral nutrients contain one or more metal ions
    (e.g., Ca, Fe, Mg, etc.). The function checks if the molecule contains any metallic
    atom from a predefined set. (This is a heuristic and does not cover every case).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is recognized as a mineral nutrient by our heuristic, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a set of metal elements commonly encountered in nutritional minerals.
    # This set is based on the examples given.
    metal_set = {"Pd", "K", "Fe", "Cs", "Ca", "Zn", "Al", "Mg", "Sb", "Ba", "Na", "La"}
    
    # Count the metal atoms from our list found in the molecule.
    found_metals = set()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in metal_set:
            found_metals.add(symbol)
    
    # If no metals from our list are found, we are unlikely to have a mineral nutrient.
    if not found_metals:
        return False, "No recognized mineral nutrient metal ions found"
    
    # Optionally, one could further require the presence of additional inorganic fragments,
    # such as oxyanions (e.g. phosphate, sulfate etc). For this simple heuristic, we assume that
    # the presence of at least one metal from the mineral_set is sufficient.
    
    metals_list = ", ".join(sorted(found_metals))
    reason = f"Found metal nutrient ion(s): {metals_list}"
    return True, reason

# Example usage (for testing):
if __name__ == "__main__":
    # Test examples from the prompt (each SMILES corresponds to a mineral nutrient)
    test_smiles = {
        "Potassium hexachloropalladate(IV)": "[Pd-2](Cl)(Cl)(Cl)(Cl)(Cl)Cl.[K+].[K+]",
        "Iron(3+) phosphate": "[Fe+3].[O-]P([O-])(=O)[O-]",
        "Caesium formate": "[Cs+].[H]C([O-])=O",
        "Calcium silicate": "[Ca++].[Ca++].[O-][Si]([O-])([O-])[O-]",
        "Zinc nitrate": "[Zn++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",
        "Aluminium sulfate octadecahydrate": "O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O",
        "Magnesium distearate": "[Mg++].CCCCCCCCCCCCCCCCCC([O-])=O.CCCCCCCCCCCCCCCCCC([O-])=O",
        "Antimony pentafluoride": "[Sb](F)(F)(F)(F)F",
        "Caesium chloride": "[Cl-].[Cs+]",
        "Calcium hypochlorite": "Cl[O-].[Ca+2].Cl[O-]",
        "Lanthanum trichloride": "Cl[La](Cl)Cl",
        "Magnesium phosphate": "P([O-])([O-])([O-])=O.P([O-])([O-])([O-])=O.[Mg+2].[Mg+2].[Mg+2]",
        "Calcium dichloride": "[Cl-].[Cl-].[Ca++]",
        "Trisodium phosphate": "[Na+].[Na+].[Na+].[O-]P([O-])([O-])=O",
        "Barium sulfate": "[O-]S([O-])(=O)=O.[Ba+2]",
        "Aluminium sulfate (anhydrous)": "[Al+3].[Al+3].[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O.[O-]S([O-])(=O)=O",
        "Barium nitrate": "[Ba++].[O-][N+]([O-])=O.[O-][N+]([O-])=O",
        "Calcium monohydroxide": "O[Ca]",
        "Calcium sulfate": "[Ca++].[O-]S([O-])(=O)=O",
        "Magnesium dipropionate": "[Mg++].CCC([O-])=O.CCC([O-])=O",
        "Tricalcium bis(phosphate)": "[Ca++].[Ca++].[Ca++].[O-]P([O-])([O-])=O.[O-]P([O-])([O-])=O",
        "Disodium hydrogenphosphate": "[Na+].[Na+].OP([O-])([O-])=O",
        "Potassium chloride": "[Cl-].[K+]",
        "Calcium dihydroxide": "[OH-].[OH-].[Ca++]",
        "Magnesium dichloride": "[Mg++].[Cl-].[Cl-]",
        "Magnesium sulfate": "[Mg++].[O-]S([O-])(=O)=O",
        "Barium acetate": "[Ba++].CC([O-])=O.CC([O-])=O",
        "Barium carbonate": "[Ba++].[O-]C([O-])=O",
        "Calcium carbonate": "[Ca+2].C(=O)([O-])[O-]",
        "Calcium difluoride": "[F-].[F-].[Ca++]",
        "Potassium sulfate": "[K+].[K+].[O-]S([O-])(=O)=O",
        "Calcium hydrogenphosphate": "[Ca++].[H]OP([O-])([O-])=O"
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_mineral_nutrient(smi)
        print(f"{name}: {result}, {reason}")