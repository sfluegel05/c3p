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

    # Initialize components list
    components = []
    
    # Check for peptide/amino acid backbone
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]")
    if mol.HasSubstructMatch(amino_acid_pattern):
        aa_matches = len(mol.GetSubstructMatches(amino_acid_pattern))
        if aa_matches >= 1:
            components.append(f"peptide/amino acid chain ({aa_matches} residues)")

    # Check for CoA core
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[CH]([OH])C(C)(C)COP([OH])(=O)OP([OH])(=O)OCC1OC([CH]C(O)C1OP([OH])([OH])=O)n1cnc2c(N)ncnc12")
    if mol.HasSubstructMatch(coa_pattern):
        components.append("CoA")

    # Check for non-CoA nucleotide
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if mol.HasSubstructMatch(nucleotide_pattern) and not mol.HasSubstructMatch(coa_pattern):
        components.append("nucleotide")

    # Check for fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]C(=O)[OH,O-,N,S]")
    if mol.HasSubstructMatch(fatty_acid_pattern):
        if rdMolDescriptors.CalcNumRotatableBonds(mol) > 8:
            components.append("fatty acid")

    # Check for glutathione core
    gsh_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H](CCC(=O)[NX3,NX4+][CX4H](CS)[CX3](=O)[NX3,NX4+][CX4H])[CX3](=O)[O-,OH]")
    if mol.HasSubstructMatch(gsh_pattern):
        components.append("glutathione")

    # Check for significant biological modifications that count as components
    mod_patterns = {
        "indole_conjugate": ("c1ccc2[nH]ccc2c1CC", "indole derivative"),
        "DOPA_conjugate": ("c1cc(O)c(O)cc1CC[NH2]", "DOPA derivative"),
        "selenium_conjugate": ("[Se]CC[CH]([NH2])C(=O)", "selenium-containing amino acid"),
    }
    
    for pattern, (smarts, desc) in mod_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            components.append(desc)

    # Special cases for known bioconjugate types
    if "CoA" in components and "fatty acid" in components:
        return True, "Bioconjugate containing: CoA conjugated with fatty acid"
    
    if "glutathione" in components and mol.HasSubstructMatch(Chem.MolFromSmarts("SC[CH]([NH2])C(=O)")):
        # Check if glutathione is conjugated to something
        if any(c not in ["glutathione", "peptide/amino acid chain"] for c in components):
            return True, "Bioconjugate containing: glutathione conjugate with " + "; ".join(c for c in components if c != "glutathione")

    # Count unique biological components
    unique_components = set(components)
    
    # Require at least 2 distinct major biological components
    major_components = [c for c in unique_components if any(x in c.lower() for x in 
                      ['peptide', 'nucleotide', 'coa', 'fatty acid', 'glutathione'])]
    
    if len(major_components) >= 2:
        return True, "Bioconjugate containing: " + "; ".join(major_components)
    
    # Handle case where only one major component is found
    if len(major_components) == 1:
        return False, f"Only one component found: {major_components[0]}"
        
    return False, "Insufficient biological components for bioconjugate classification"