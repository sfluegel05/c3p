"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import ChiralType

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    L-amino acids have S configuration at the alpha carbon (except cysteine and derivatives).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate 3D coordinates for chirality checking
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        pass  # Continue even if 3D generation fails

    # SMARTS patterns for alpha-amino acid core
    patterns = [
        # Various forms of alpha-amino acids including charged states and substituted amines
        '[NX3,NX4+][CX4;H1,H2]([*])[CX3](=[OX1])[OX1-,OX2H1]',
        '[NX3,NX4+][CX4;H1,H2]([*])[CX3](=[OX1])[OX2][*]',  # Ester or other derivatives
        '[NX3,NX4+;H0,H1,H2,H3][CX4;H1,H2]([*])[CX3](=[OX1])[OX1-,OX2H1,OX2]'  # More general pattern
    ]

    for pattern in patterns:
        aa_pattern = Chem.MolFromSmarts(pattern)
        if aa_pattern is None:
            continue
            
        matches = mol.GetSubstructMatches(aa_pattern)
        for match in matches:
            # Get the alpha carbon (second atom in pattern)
            alpha_carbon = mol.GetAtomWithIdx(match[1])
            
            # Must be chiral
            if alpha_carbon.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                continue
                
            # Get neighbors
            neighbors = list(alpha_carbon.GetNeighbors())
            if len(neighbors) != 4:
                continue
                
            # Classify neighbor types
            n_types = {'N': None, 'C': None, 'H': None, 'R': None}
            for n in neighbors:
                if n.GetAtomicNum() == 7:  # N
                    n_types['N'] = n
                elif n.GetAtomicNum() == 1:  # H
                    n_types['H'] = n
                elif n.GetAtomicNum() == 6 and any(a.GetAtomicNum() == 8 for a in n.GetNeighbors()):  # C(=O)
                    n_types['C'] = n
                else:
                    n_types['R'] = n
                    
            # Must have all required groups
            if not all(n_types[k] is not None for k in ['N', 'C', 'H']):
                continue

            # Check configuration
            # For L-amino acids:
            # - Looking from H: N, R, COOH should be counterclockwise (S configuration)
            # - Exception: If R group has higher CIP priority than COOH, should be clockwise (R configuration)
            
            # Use RDKit's assignStereochemistry to get the correct configuration
            Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
            
            # L-amino acids typically have S configuration (except cysteine derivatives)
            # Check for sulfur-containing R group which would invert the expected configuration
            has_sulfur = False
            if n_types['R'] is not None:
                for atom in mol.GetAtoms():
                    if atom.GetAtomicNum() == 16:  # Sulfur
                        has_sulfur = True
                        break
            
            expected_tag = '@' if has_sulfur else '@@'
            atom_smiles = Chem.MolToSmiles(mol, rootedAtAtom=match[1])
            
            if (has_sulfur and '@' in atom_smiles and '@@' not in atom_smiles) or \
               (not has_sulfur and '@@' in atom_smiles):
                return True, "Found L-alpha-amino acid with correct stereochemistry"
            
    # If we get here, check if we found any amino acid pattern
    for pattern in patterns:
        aa_pattern = Chem.MolFromSmarts(pattern)
        if aa_pattern and mol.HasSubstructMatch(aa_pattern):
            return False, "Found alpha-amino acid pattern but incorrect or unspecified stereochemistry"
    
    return False, "No alpha-amino acid pattern found"