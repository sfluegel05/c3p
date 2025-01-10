"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Classical cannabinoid patterns (more general)
    thc_core = Chem.MolFromSmarts("c1c(O)cc(CCCCC)c2c1OC(C)(C)C1CCC(C)=CC21") # THC core
    cbd_core = Chem.MolFromSmarts("c1c(O)cc(CCCCC)c(O)c1CC1C=C(C)CCC1") # CBD core
    cbg_core = Chem.MolFromSmarts("c1c(O)cc(CCCCC)c(O)c1CC=C(C)CCC=C(C)C") # CBG core
    
    # Cannabinoid acid patterns
    thca_core = Chem.MolFromSmarts("c1c(O)c(C(=O)O)c(CCCCC)c2c1OC(C)(C)C1CCC(C)=CC21")
    cbda_core = Chem.MolFromSmarts("c1c(O)c(C(=O)O)c(CCCCC)c(O)c1CC1C=C(C)CCC1")
    cbga_core = Chem.MolFromSmarts("c1c(O)c(C(=O)O)c(CCCCC)c(O)c1CC=C(C)CCC=C(C)C")
    
    # Synthetic cannabinoid patterns
    indole_core = Chem.MolFromSmarts("c1ccc2c(c1)c(C(=O)[#6])n([CH2][CH2][CH2][CH2][#6])2") # More specific indole
    cp_core = Chem.MolFromSmarts("c1c(O)c([CH2,CH]C2CC[CH](O)CC2)cc(C(C)(C)[CH2,CH3])c1") # CP-like
    
    # Endocannabinoid patterns
    arachidonic_chain = Chem.MolFromSmarts("CCCCC=CC=CC=CC=CC=CC") # Arachidonic acid chain
    glycerol_ester = Chem.MolFromSmarts("OCC(O)COC(=O)[CH2][CH2]") # Glycerol ester
    ethanolamine = Chem.MolFromSmarts("C(=O)NCCO") # Ethanolamine

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Count oxygens (cannabinoids typically have 2-3 oxygens)
    o_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[O]")))
    
    # Classical cannabinoids and their acids
    classical_cores = [
        (thc_core, "THC-like"), 
        (cbd_core, "CBD-like"),
        (cbg_core, "CBG-like"),
        (thca_core, "THCA-like"),
        (cbda_core, "CBDA-like"),
        (cbga_core, "CBGA-like")
    ]
    
    for core, core_type in classical_cores:
        if mol.HasSubstructMatch(core):
            if 280 < mol_wt < 400 and ring_count >= 2 and 2 <= o_count <= 4:
                return True, f"Classical cannabinoid with {core_type} core structure"

    # Synthetic cannabinoids
    if mol.HasSubstructMatch(indole_core):
        if 300 < mol_wt < 450 and ring_count >= 3:
            return True, "Synthetic cannabinoid with indole core"
            
    if mol.HasSubstructMatch(cp_core):
        if 300 < mol_wt < 450 and ring_count >= 2:
            return True, "Synthetic cannabinoid with CP-like core"

    # Endocannabinoids
    has_arachidonic = mol.HasSubstructMatch(arachidonic_chain)
    has_glycerol = mol.HasSubstructMatch(glycerol_ester)
    has_ethanolamine = mol.HasSubstructMatch(ethanolamine)
    
    if has_arachidonic:
        if has_glycerol:
            if 300 < mol_wt < 500 and rotatable_bonds > 10 and o_count == 4:
                return True, "2-arachidonoyl glycerol (2-AG) type endocannabinoid"
        if has_ethanolamine:
            if 300 < mol_wt < 500 and rotatable_bonds > 10 and o_count == 3:
                return True, "Anandamide-type endocannabinoid"

    return False, "Does not match cannabinoid structural patterns"