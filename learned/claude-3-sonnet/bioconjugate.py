"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: bioconjugate
A molecular entity consisting of at least 2 biological molecules covalently linked together.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate must contain at least 2 biological molecules covalently linked.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize components list and biological modification flags
    components = []
    has_biological_modification = False
    
    # Check for peptide/amino acid backbone
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]")
    if mol.HasSubstructMatch(amino_acid_pattern):
        aa_matches = len(mol.GetSubstructMatches(amino_acid_pattern))
        if aa_matches >= 1:
            components.append(f"peptide/amino acid chain ({aa_matches} residues)")

    # Check for nucleotide core (excluding CoA)
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if mol.HasSubstructMatch(nucleotide_pattern):
        # Only count if not part of CoA
        coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[CH]([OH])C(C)(C)COP([OH])(=O)OP([OH])(=O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP([OH])([OH])=O")
        if not mol.HasSubstructMatch(coa_pattern):
            components.append("nucleotide")

    # Check for significant biological modifications
    mod_patterns = {
        "indole": ("c1ccc2[nH]ccc2c1", "indole modification"),
        "IAN": ("C(#N)Cc1c[nH]c2ccccc12", "IAN group"),
        "DOPA": ("c1cc(O)c(O)cc1CC", "DOPA modification"),
        "selenium": ("[Se]", "selenium modification"),
        "thiol_conjugate": ("SC[CH]([NH2])C(=O)", "thiol conjugate")
    }
    
    for pattern, (smarts, desc) in mod_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            has_biological_modification = True
            components.append(desc)

    # Check for lipid/fatty acid conjugation
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]C(=O)[OH,O-,N]")
    if mol.HasSubstructMatch(fatty_acid_pattern):
        # Verify chain length to avoid small modifications
        if rdMolDescriptors.CalcNumRotatableBonds(mol) > 8:
            components.append("fatty acid")

    # Special handling for glutathione conjugates
    gsh_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H](CCC(=O)[NX3,NX4+][CX4H](CS)[CX3](=O)[NX3,NX4+][CX4H])[CX3](=O)[O-,OH]")
    if mol.HasSubstructMatch(gsh_pattern):
        # Look for conjugation to the cysteine thiol
        if mol.HasSubstructMatch(Chem.MolFromSmarts("SCC[CH]([NH2])C(=O)")):
            has_biological_modification = True

    # Final decision logic
    unique_components = set(components)
    
    # Case 1: Multiple distinct biological components
    if len(unique_components) >= 2:
        return True, "Bioconjugate containing: " + "; ".join(unique_components)
    
    # Case 2: One component with significant biological modification
    elif len(unique_components) == 1 and has_biological_modification:
        return True, f"Bioconjugate containing: {list(unique_components)[0]} with biological modification"
    
    # Case 3: Single component without significant modification
    elif len(unique_components) == 1:
        return False, f"Only one component found: {list(unique_components)[0]}"
    
    # Case 4: No biological components
    else:
        return False, "No biological components identified"