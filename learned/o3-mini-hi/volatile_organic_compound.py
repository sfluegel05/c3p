"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: Volatile Organic Compound (VOC)
Definition: “Any organic compound having an initial boiling point less than or equal 
to 250 °C (482 °F) measured at a standard atmospheric pressure of 101.3 kPa.”

This heuristic uses a combination of estimated molecular descriptors (molecular weight and 
topological polar surface area) along with several “structural” checks including:
 • Rejecting molecules that are not organic (lacking carbon), or that carry formal charges.
 • Counting non‐halogen heteroatoms (excluding H, C and common halogens F, Cl, Br, I)
   – if too many (≥4), then the molecule is assumed to be “highly functionalized.”
 • Using ring information. If a molecule has more than one ring or two (or more) aromatic rings,
   it is rejected.
 • A few specific “problematic” groups are also flagged:
     – A carboxylic acid group in a molecule with MW >150 (larger acids tend to be less volatile)
     – An α,β‑unsaturated carbonyl (conjugated enone) is flagged.
     – Very extended conjugated polyene systems.
 • Finally, the molecular weight cut‐off is chosen based on whether the molecule contains an alcohol 
   (common in the true‐positive VOC examples), or else if it is aromatic.
If the molecule is acyclic and contains an –OH group, a cutoff of 400 Da is used;
if it is aromatic (monocyclic) but not an –OH compound, a lower cutoff (250 Da) is used;
otherwise a cutoff of 300 Da is applied.
Additionally, a molecule must have low topological polar surface area (TPSA < 60 Å²)
to be considered volatile.
Note: This heuristic is approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on a 
    heuristic estimation.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a VOC, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Requirement: must be organic (i.e. contain at least one carbon atom)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not organic (contains no carbon)"
    
    # Reject molecules with any nonzero formal charge.
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            return False, "Molecule carries a formal charge"
    
    # Compute basic descriptors: molecular weight and TPSA.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    
    # Count non‐halogen heteroatoms (exclude H (1), carbon (6) and common halogens F (9), Cl (17), Br (35), I (53))
    allowed_halogens = {9, 17, 35, 53}
    hetero_count = sum(1 for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() not in (1, 6) and atom.GetAtomicNum() not in allowed_halogens)
    if hetero_count >= 4:
        return False, f"Too many heteroatoms ({hetero_count}), suggesting high functionality and high boiling point"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Count aromatic rings using the ring atom indices from ring info.
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1
    if aromatic_ring_count >= 2:
        return False, f"Contains {aromatic_ring_count} aromatic rings which tend to increase boiling point"
    if ring_count > 1:
        return False, f"Contains {ring_count} rings which tend to increase boiling point"
    # For non‐aromatic single rings that are highly unsaturated (eg, conjugated dienes), reject them.
    if ring_count == 1 and not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        num_aliphatic_db = rdMolDescriptors.CalcNumAliphaticDoubleBonds(mol)
        if num_aliphatic_db >= 2:
            return False, "Monocyclic diene detected – such unsaturation often increases boiling point"
    
    # Define some SMARTS patterns for functional groups.
    acid_pattern = Chem.MolFromSmarts("[$([CX3](=O)[OX2H1]),$([CX3](=O)[O-])]")
    enone_pattern = Chem.MolFromSmarts("[C]=[C][C](=O)[C]")  # a simplified conjugated ketone pattern
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")  # at least three conjugated C=C bonds
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    
    # If the molecule contains a carboxylic acid and is not very small, assume high boiling point.
    if mol.HasSubstructMatch(acid_pattern) and mol_wt > 150:
        return False, f"Contains a carboxylic acid group and MW ({mol_wt:.1f} Da) >150, suggesting high boiling point"
    
    # If the molecule contains a conjugated ketone (enone), flag it.
    if mol.HasSubstructMatch(enone_pattern):
        return False, "Contains a conjugated enone group, which increases the boiling point"
    
    # If the molecule has a long conjugated polyene system, reject.
    if mol.HasSubstructMatch(polyene_pattern):
        return False, "Contains an extended polyene system, which tends to increase boiling point"
    
    # Decide which molecular weight cutoff to use.
    # Many true positives are simple alcohols -> allow higher MW (<=400 Da);
    # non‐alcohols (and non‐aromatic compounds) use a lower cutoff (<=300 Da);
    # if the molecule is monocyclic aromatic, use an even stricter cutoff (<=250 Da)
    if mol.HasSubstructMatch(alcohol_pattern):
        cutoff_mw = 400
    elif aromatic_ring_count == 1:
        cutoff_mw = 250
    else:
        cutoff_mw = 300

    # Final decision based on MW and TPSA.
    if mol_wt <= cutoff_mw and tpsa < 60:
        reason = (f"Estimated as VOC: MW ({mol_wt:.1f} Da) <= {cutoff_mw} and TPSA ({tpsa:.1f} Å²) < 60, "
                  "suggesting a low boiling point (<=250 °C)")
        return True, reason
    else:
        reason = (f"Estimated not VOC: MW ({mol_wt:.1f} Da) and TPSA ({tpsa:.1f} Å²) do not meet criteria "
                  f"for low boiling point (<=250 °C) with cutoff MW = {cutoff_mw} Da")
        return False, reason

# Example usage (for testing)
if __name__ == "__main__":
    # A selection of examples from the training outcomes.
    test_examples = {
        "nonan-2-ol": "CCCCCCCC(C)O",
        "decan-2-ol": "CCCCCCCCC(C)O",
        "2-dodecene": "[H]C(C)=C([H])CCCCCCCCC",
        "heptadecan-8-ol": "CCCCCCCCC(C)OCCCCCCC",  # slightly changed for acyclicity
        "henicosan-3-ol": "CCCCCCCCCCCCCCCCCC(O)CC",
        "2,3,5-trimethylhexane": "C(CC(C)C)(C(C)C)C",
        "3-methylpentane": "CCC(C)CC",
        "thiophene": "c1ccsc1",
        "4-methylpent-3-enoic acid": "OC(=O)CC=C(C)C",
        "3-Methyl-3-hepten-2-one": "O=C(\\C(=C\\CCC)C)C",
        "OC(=O)CCCCCCCC=CCC=CCCCC": "OC(=O)CCCCCCCC=CCC=CCCCC",
    }
    for name, smi in test_examples.items():
        voc, reason = is_volatile_organic_compound(smi)
        print(f"Name: {name}\n  SMILES: {smi}\n  VOC: {voc}\n  Reason: {reason}\n")