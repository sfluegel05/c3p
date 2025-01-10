"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Convert SMILES to RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expanded SMARTS patterns for purine-like and pyrimidine-like structures
    purine_like_smarts = [
        Chem.MolFromSmarts("c1ncnc2[nH]cnn12"),  # Basic purine extended with variations
        Chem.MolFromSmarts("c1nc[nH]c2c(ncnc12)"),  # Substituted purine ring
    ]
    
    pyrimidine_like_smarts = [
        Chem.MolFromSmarts("c1[nH]cncnc1"),       # Basic pyrimidine with variations
        Chem.MolFromSmarts("c1c[nH]cnc[nH]c1"),   # Alterations with additional rings
        Chem.MolFromSmarts("c1ncc(=O)[nH]n1"),    # Typical keto and imidazole derivatives
        Chem.MolFromSmarts("c1nc[nH]cnc1"),       # Expanded pyrimidine scaffold
        Chem.MolFromSmarts("c1[nH]cc(=O)[nH]c1"), # Variants with keto positions
    ]
    
    # Check for purine-like or pyrimidine-like structures
    for pattern in purine_like_smarts:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches purine-like structure"

    for pattern in pyrimidine_like_smarts:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches pyrimidine-like structure"
    
    # Check for functional groups characteristic of nucleobase analogues
    functional_groups = [
        Chem.MolFromSmarts("[CX3]=[OX1]"),    # Keto/Carbonyl group
        Chem.MolFromSmarts("[NX3;H2]"),       # Amino group
        Chem.MolFromSmarts("[OX2H]"),         # Hydroxy group
        Chem.MolFromSmarts("[#16]"),          # Incorporation of sulfur (thio)
        Chem.MolFromSmarts("[Cl,Br,I,F]"),    # Halogens
    ]
    
    # Consider complexity: require at least one functional group and a basic ring
    if any(mol.HasSubstructMatch(fg) for fg in functional_groups):
        if mol.GetNumAtoms() > 6:  # Filter small molecules
            if metadata_for_atom_rings(mol):
                return True, "Contains significant nucleobase-related functional groups and structure"
    
    # None matched
    return False, "No nucleobase analogue characteristics found"

def metadata_for_atom_rings(mol):
    """
    Check molecule for metadata indicative of nucleobases.
    For example, 3-atom cons mostly Nitrogen and Carbon.
    """
    ring_info = mol.GetRingInfo()
    for atom_counts in ring_info.AtomRings():
        if len(atom_counts) <= 6:
            # Check number of Nitrogens in small rings
            n_count = sum(1 for i in atom_counts if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
            if n_count >= 2:  # Purines and related structures often have 2+
                return True
    return False