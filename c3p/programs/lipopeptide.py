"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a lipopeptide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of long aliphatic chain as part of lipid moiety
    # Let's check for paths of length 12 to 24, which are representative of typical fatty acids
    long_chain_found = False
    for length in range(12, 25):
        num_chains = len(rdmolops.FindAllPathsOfLengthN(mol, length, useBonds=True))
        if num_chains > 0:
            long_chain_found = True
            break

    if not long_chain_found:
        return False, "No long hydrocarbon chain found (indicating lipid motif)"

    # Check for presence of amide bonds typically found in peptide chains
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 3:  # Increased the requirement to a minimum of 3
        return False, f"Found {len(amide_matches)} amide bonds, need at least 3 for peptide linkage"

    # Count chiral centers as a proxy for peptide sequence complexity
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 2:  # Reduced to 2, as lipopeptides may be simple in structure
        return False, f"Found {len(chiral_centers)} chiral centers, indicating insufficient peptide component"

    # This additional pattern can help identify more peptide-like properties
    # e.g., looking for typical N-terminus
    # nitrogen_pattern = Chem.MolFromSmarts("[NH2,NH3+,NH+]")

    return True, "Contains long hydrocarbon chain and multiple amide bonds and chiral centers indicative of lipopeptide structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46895',
                          'name': 'lipopeptide',
                          'definition': 'A compound consisting of a peptide '
                                        'with attached lipid.',
                          'parents': ['CHEBI:16670', 'CHEBI:18059'],
                          'xrefs': [   'PMID:19889045',
                                       'PMID:20545290',
                                       'PMID:23131643',
                                       'PMID:23318669',
                                       'Wikipedia:Lipopeptide'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.5302325581395348 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               '[C@@]12(N(C[C@@H](C1)O)C([C@H]([C@@H](C)O)NC(=O)[C@](C[C@H]([C@@H](O)NC([C@@H]3[C@H](CCN3C([C@H]([C@@H](CC(N)=O)O)NC(=O)[C@H]([C@@H]([C@H](C4=CC=C(C=C4)O)O)O)NC2=O)=O)O)=O)O)(NC(CCCCCCCC[C@H](C[C@H](CC)C)C)=O)[H])=O)[H] '
               'NAME: pneumocandin B0 REASON: CORRECT Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)CC(=O)N[C@H](Cc1ccccc1)C(=O)N[C@H]([C@@H](C)OC1O[C@@H](C)[C@@H](O)[C@@H](O)[C@H]1O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O)C(=O)N[C@H](C)C(=O)N[C@@H](C)COC1O[C@@H](C)[C@H](OC)[C@@H](OC)[C@H]1O '
               'NAME: GPL-1 REASON: CORRECT Contains long hydrocarbon chain '
               'and multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'CC(C)[C@H](N(C)C(=O)[C@H](C(C)C)N(C)C(=O)[C@H](C(C)C)N(C)C(=O)[C@@H](C)CCCCC#C)C(=O)N(C)[C@@H](Cc1ccccc1)C(N)=O '
               'NAME: dragonamide A REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)CCc1ccccc1C(F)(F)F)NC(=O)[C@H](CC(O)=O)NC(C)=O)[C@@H](C)O)C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O '
               'NAME: '
               'Ac-Asp-N(6)-{3-[2-(trifluoromethyl)phenyl]propanoyl}-KATIGFEVQEE '
               'REASON: CORRECT Contains long hydrocarbon chain and multiple '
               'amide bonds and chiral centers indicative of lipopeptide '
               'structure\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCC(CC(=O)N[C@H](Cc1ccccc1)C(=O)N[C@H]([C@@H](C)O[C@@H]1O[C@@H](C)[C@@H](OC(C)=O)[C@@H](OC(C)=O)[C@H]1O)C(=O)N[C@H](C)C(=O)N[C@@H](C)CO[C@@H]1O[C@@H](C)[C@H](OC)[C@@H](OC)[C@H]1O)OC '
               'NAME: GPL C-II REASON: CORRECT Contains long hydrocarbon chain '
               'and multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)c1cccc(c1)C(F)(F)F)NC(=O)[C@H](CC(O)=O)NC(C)=O)[C@@H](C)O)C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O '
               'NAME: Ac-Asp-N(6)-[3-(trifluoromethyl)benzoyl]-KATIGFEVQEE '
               'REASON: CORRECT Contains long hydrocarbon chain and multiple '
               'amide bonds and chiral centers indicative of lipopeptide '
               'structure\n'
               ' * SMILES: '
               'CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)\\C=C\\c1cccc(c1)[N+]([O-])=O)NC(=O)[C@H](CC(O)=O)NC(C)=O)[C@@H](C)O)C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O '
               'NAME: '
               'Ac-Asp-N(6)-[(2E)-3-(3-nitrophenyl)prop-2-enoyl]-KATIGFEVQEE '
               'REASON: CORRECT Contains long hydrocarbon chain and multiple '
               'amide bonds and chiral centers indicative of lipopeptide '
               'structure\n'
               ' * SMILES: '
               '[C@@]12(N(C[C@@H](C1)O)C([C@H]([C@@H](C)O)NC(=O)[C@](C[C@H]([C@@H](NCCN)NC([C@@H]3[C@H](CCN3C([C@H]([C@@H](CCN)O)NC(=O)[C@H]([C@@H]([C@H](C4=CC=C(C=C4)O)O)O)NC2=O)=O)O)=O)O)(NC(CCCCCCCC[C@H](C[C@H](CC)C)C)=O)[H])=O)[H] '
               'NAME: caspofungin REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'O.CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)CC(=O)N[C@H](Cc1ccccc1)C(=O)N[C@H]([C@@H](C)OC1O[C@@H](C)[C@@H](O)[C@@H](O)[C@H]1O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O[C@@H]2O[C@@H]3CO[C@@](C)(O[C@H]3[C@H](OC)[C@H]2O)C(O)=O)[C@H]1O)C(=O)N[C@H](C)C(=O)N[C@@H](C)COC1O[C@@H](C)[C@H](OC)[C@@H](OC)[C@H]1O '
               'NAME: GPL-21 REASON: CORRECT Contains long hydrocarbon chain '
               'and multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)NCCCC[C@@H](NC(=O)CC[C@@H](NC(=O)[C@H](C)NC(=O)[C@@H](C)O[C@@H]1[C@@H](NC(C)=O)[C@H](O)O[C@H](CO)[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O)C(N)=O)C(O)=O '
               'NAME: beta-D-GlcNAc-(1->4)-MDP-Lys(L18) REASON: CORRECT '
               'Contains long hydrocarbon chain and multiple amide bonds and '
               'chiral centers indicative of lipopeptide structure\n'
               ' * SMILES: '
               'CCCCCCCCOc1ccc(cc1)C(=O)N[C@H]1C[C@@H](O)[C@@H](O)NC(=O)[C@@H]2[C@@H](O)[C@@H](C)CN2C(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC1=O)[C@@H](C)O)[C@H](O)[C@@H](O)c1ccc(O)cc1)[C@@H](C)O '
               'NAME: cilofungin REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'C(N[C@@H]1C(N[C@@H](C(N[C@@H](CCCCN)C(N[C@H](C(N[C@H](C(NCC(N[C@@H](C(N[C@H](C(N2[C@H](C(NC1C)=O)C[C@H](C2)C)=O)C(C)C)=O)C(C)C(=O)O)=O)=O)CC(O)=O)=O)C(C(O)=O)O)=O)=O)C(C)C)=O)([C@@H](NC(/C=C/C=C\\CCC(C)C)=O)C(C)C(=O)O)=O '
               'NAME: malacidin A REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'CCCCCCCC(=O)NCCCC[C@H](NC(=O)[C@H](CC(O)=O)NC(C)=O)C(=O)N[C@@H](C)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H]([C@@H](C)CC)C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O '
               'NAME: Ac-Asp-N(6)-(octanoyl)-KATIGFEVQEE REASON: CORRECT '
               'Contains long hydrocarbon chain and multiple amide bonds and '
               'chiral centers indicative of lipopeptide structure\n'
               ' * SMILES: '
               'CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)\\C=C\\c1cccc(c1)C(F)(F)F)NC(=O)[C@H](CC(O)=O)NC(C)=O)[C@@H](C)O)C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O '
               'NAME: '
               'Ac-Asp-N(6)-{(2E)-3-[3-(trifluoromethyl)phenyl]prop-2-enoyl}-KATIGFEVQEE '
               'REASON: CORRECT Contains long hydrocarbon chain and multiple '
               'amide bonds and chiral centers indicative of lipopeptide '
               'structure\n'
               ' * SMILES: '
               '[H][C@@]1(CCCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1 '
               'NAME: surfactin D REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'CC[C@@H](C)CCCCC(=O)N[C@@H](CCN)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CCN)C(=O)N[C@H]1CCNC(=O)[C@@H](NC(=O)[C@H](CCN)NC(=O)[C@H](CCN)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](CC(C)C)NC(=O)[C@H](CCN)NC1=O)[C@@H](C)O '
               'NAME: colistin A REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'CC[C@H](C)[C@H](NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CCCCNC(=O)COCCOCCNC(=O)COCCOCCNC(=O)CC[C@@H](NC(=O)CCCCCCCCCCCCCCCCC(O)=O)C(O)=O)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CCC(N)=O)NC(=O)CNC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC1=CC=C(O)C=C1)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(O)=O)NC(=O)C(C)(C)NC(=O)[C@@H](N)CC1=CN=CN1)[C@@H](C)O)[C@@H](C)O)C(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](CC1=CNC2=C1C=CC=C2)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCCNC(N)=N)C(=O)NCC(=O)N[C@@H](CCCNC(N)=N)C(=O)NCC(O)=O '
               'NAME: semaglutide REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'CCCCCC[C@H](C)[C@@H]1CC(=O)NCC(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O1 '
               'NAME: emericellamide D REASON: CORRECT Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               '[H][C@@]1(CCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O1 '
               'NAME: surfactin A REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)\\C=C\\c1ccc(Cl)cc1)NC(=O)[C@H](CC(O)=O)NC(C)=O)[C@@H](C)O)C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O '
               'NAME: '
               'Ac-Asp-N(6)-[(2E)-3-(4-chlorophenyl)prop-2-enoyl]-KATIGFEVQEE '
               'REASON: CORRECT Contains long hydrocarbon chain and multiple '
               'amide bonds and chiral centers indicative of lipopeptide '
               'structure\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)CC(=O)N[C@H](Cc1ccccc1)C(=O)N[C@H]([C@@H](C)OC1O[C@@H](C)[C@@H](O)[C@@H](O)[C@H]1O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O[C@@H]2O[C@@H](C)[C@@H](OC(C)=O)[C@@H](OC)[C@@H]2OC)[C@H]1O)C(=O)N[C@H](C)C(=O)N[C@@H](C)COC1O[C@@H](C)[C@H](OC)[C@@H](OC)[C@H]1O '
               'NAME: acetyl-GPL-2 REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)CC(=O)N[C@H](Cc1ccccc1)C(=O)N[C@H]([C@@H](C)OC1O[C@@H](C)[C@@H](O)[C@@H](O)[C@H]1O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O[C@@H]2O[C@@H](C)[C@@H](O[C@@H]3O[C@@H]([C@@H](O[C@@H]4O[C@@H](C)[C@@H](OC(C)=O)[C@@H](OC)[C@@H]4OC)[C@H](O)[C@H]3O)C(O)=O)[C@@H](OC)[C@@H]2OC)[C@H]1O)C(=O)N[C@H](C)C(=O)N[C@@H](C)COC1O[C@@H](C)[C@H](OC)[C@@H](OC)[C@H]1O '
               'NAME: GPL-9I REASON: CORRECT Contains long hydrocarbon chain '
               'and multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'CCCCCC[C@H](C)[C@H]1OC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)CNC(=O)[C@@H]1C)C(C)C '
               'NAME: emericellamide A REASON: CORRECT Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'CCCCCCCCOc1ccc(cc1)-c1ccc(cc1)C(=O)N[C@H]1CC(CNC(=O)[C@@H]2[C@@H](O)[C@@H](C)CN2C(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC1=O)[C@@H](C)O)[C@H](O)Cc1ccc(O)cc1)NCCN '
               'NAME: aminocandin REASON: CORRECT Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'C(CCCCCCCC)CCCCC(N[C@H](C(=O)N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(=O)N[C@H](C(=O)O)CCCCN)CCCCN)=O)CCCCN)=O)CCCCN)=O)CO)CSCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC)=O '
               'NAME: '
               'N-palmitoyl-S-[2,3-bis(palmitoyloxy)propyl]-Cys-Ser-Lys-Lys-Lys-Lys '
               'REASON: CORRECT Contains long hydrocarbon chain and multiple '
               'amide bonds and chiral centers indicative of lipopeptide '
               'structure\n'
               'False positives: SMILES: '
               'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C '
               'NAME: Cyanopeptolin D REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C '
               'NAME: Ergovaline REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O '
               'NAME: desferrialbomycin epsilon(3-) REASON: WRONGLY CLASSIFIED '
               'Contains long hydrocarbon chain and multiple amide bonds and '
               'chiral centers indicative of lipopeptide structure\n'
               ' * SMILES: '
               'C[C@@H](C1=CC=CC=C1)NC(=O)C[C@@H]2CC[C@H]([C@H](O2)CO)NC(=O)C3CCC3 '
               'NAME: '
               'N-[(2S,3R,6S)-2-(hydroxymethyl)-6-[2-oxo-2-[[(1S)-1-phenylethyl]amino]ethyl]-3-oxanyl]cyclobutanecarboxamide '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@H]3O[C@@H]([C@@H](O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4NC(=O)C)CO)[C@H](O)[C@@H]3O)CO)[C@H](O)[C@@H]1O[C@H]5[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]5CO)O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O)CO[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O)[C@H](O)[C@H]8NC(=O)C)CO)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4R,5R,6R)-5-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-[(2R,3S,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3,4-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-3-hydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'C[C@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)C3CC3)C(=O)N(C[C@@H]1OC)C)C)C(=O)C4CCCCC4 '
               'NAME: '
               'N-[(5R,6S,9R)-8-[cyclohexyl(oxo)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]cyclopropanecarboxamide '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'C1[C@@H]2[C@H]([C@H](N2C(=O)CN1C(=O)C3=NC=CN=C3)CO)C4=CC=C(C=C4)C5=CC(=CC=C5)F '
               'NAME: '
               '(6S,7R,8S)-7-[4-(3-fluorophenyl)phenyl]-8-(hydroxymethyl)-4-[oxo(2-pyrazinyl)methyl]-1,4-diazabicyclo[4.2.0]octan-2-one '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               '[NH2+]([C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(=O)N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(=O)N[C@H](C([O-])=C)CC(C)C)C(C)C)=O)CCCC[NH3+])=O)CC(C)C)=O)C(C)C)=O)CCC(N)=O)=O)CCCC[NH3+])=O)CCC(=O)N)=O)CCCC[NH3+])=O)C(C)O)=O)CCCC[NH3+])=O)CCCC[NH3+])=O)[C@H](CC)C)=O)C(C)C)=O)CCCNC(N)=[NH2+])=O)C(C)C)=O)CCCC[NH3+])=O)CC(C)C)=O)CO)=O)CCCC[NH3+])=O)C)=O)C)CCCC[NH3+])=O)C(C)O)=O)CC(C)C)=O)CC(=O)N)=O)CO)=O)C)C(=O)CCCC[C@@H]1[C@]2([C@@](CS1)(NC(N2)=O)[H])[H] '
               'NAME: HA Binder_2nd de novo REASON: WRONGLY CLASSIFIED '
               'Contains long hydrocarbon chain and multiple amide bonds and '
               'chiral centers indicative of lipopeptide structure\n'
               ' * SMILES: '
               'O([C@@H]1[C@@H](NC(=O)C)[C@@H](O[C@@H]([C@@H]1O)CO)O[C@@H]2[C@@H](O)[C@H](O[C@@H]([C@@H]2O)CO)O[C@@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]4[C@H](O)[C@@H](O)[C@@H](O[C@@H]4CO)O)[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O '
               'NAME: '
               '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3S,4S,5R,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-6-[(2S,3R,4R,5R,6R)-3-acetamido-2-[(2R,3R,4S,5S,6R)-2-[(2R,3R,4R,5R,6S)-4,5-dihydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-5-hydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-3,5-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: WRONGLY CLASSIFIED Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'C(=O)([C@@H](N)CC1=CNC2=C1C=CC=C2)N[C@H](C(=O)N[C@H](C(=O)O)C(C)C)CC3=CNC4=C3C=CC=C4 '
               'NAME: Trp-Trp-Val REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)CC2=C(C=CC(=C2)NC(=O)CC3=CC=CC=C3)O[C@H]1CN(C)CC4=CC=C(C=C4)C(=O)O)[C@H](C)CO '
               'NAME: '
               '4-[[[(2R,3S)-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-9-[(1-oxo-2-phenylethyl)amino]-2,3,4,7-tetrahydro-1,5-benzoxazonin-2-yl]methyl-methylamino]methyl]benzoic '
               'acid REASON: WRONGLY CLASSIFIED Contains long hydrocarbon '
               'chain and multiple amide bonds and chiral centers indicative '
               'of lipopeptide structure\n'
               ' * SMILES: '
               'C1([C@@H]([C@H]1C=C(C)C)C(=O)SCCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]2O[C@@H](N3C4=C(C(=NC=N4)N)N=C3)[C@@H]([C@@H]2OP([O-])([O-])=O)O)(=O)[O-])(=O)[O-])(C)C)O)=O)(C)C '
               'NAME: (R,R)-chrysanthemoyl-CoA(4-) REASON: WRONGLY CLASSIFIED '
               'Contains long hydrocarbon chain and multiple amide bonds and '
               'chiral centers indicative of lipopeptide structure\n'
               ' * SMILES: '
               'C1C[C@H]([C@@H](O[C@@H]1CC(=O)NCCC2=CC=NC=C2)CO)NC(=O)CN3CCOCC3 '
               'NAME: '
               '2-[(2S,5R,6R)-6-(hydroxymethyl)-5-[[2-(4-morpholinyl)-1-oxoethyl]amino]-2-oxanyl]-N-(2-pyridin-4-ylethyl)acetamide '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'ClC1=C(OC)C=CC(=C1)C[C@H]2NC(=O)C=CC[C@H](OC(=O)[C@@H](OC([C@@H](CNC2=O)C)=O)CC(C)C)C[C@H]3O[C@@H]3C4=CC=CC=C4 '
               'NAME: Cryptophycin-40 REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'O=C(N1[C@H](CCC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC2=CC=CC=C2)CC3=CC=CC=C3 '
               'NAME: Phe-Phe-Pro REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'CCC(=O)NC1=CC2=C(C=C1)C(=O)N(C[C@@H]([C@@H](CN([C@@H](CO2)C)C(=O)C3CC3)C)OC)C '
               'NAME: '
               'N-[(5R,6R,9R)-8-[cyclopropyl(oxo)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]propanamide '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'O[C@@H]([C@H](NC(=O)[C@@H](N)CC1=CC=CC=C1)C(=O)N[C@@H](CC=2C=3C(NC2)=CC=CC3)C(O)=O)C '
               'NAME: Phe-Thr-Trp REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'C[C@H]1CN(C(=O)C2=C(C(=CC=C2)NC(=O)C3CCOCC3)O[C@@H]1CN(C)CC4=CC=NC=C4)[C@H](C)CO '
               'NAME: '
               'N-[(2S,3S)-5-[(2R)-1-hydroxypropan-2-yl]-3-methyl-2-[[methyl(pyridin-4-ylmethyl)amino]methyl]-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-10-yl]-4-oxanecarboxamide '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'O([C@H]1[C@@H](O[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H](O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O)[C@@H]1O)CO)[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(=O)NC '
               'NAME: '
               '(2S,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5R,6S)-3-[(2S,3R,4R,5R,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-5-hydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-4-yl]oxy-4-hydroxy-N-methyl-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxamide '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'C1COCCC1C(=O)NC2=CC3=C(C=C2)O[C@@H]4[C@H]3C[C@H](O[C@@H]4CO)CC(=O)NCC5=CN=CC=C5 '
               'NAME: '
               'N-[(1R,3S,4aS,9aR)-1-(hydroxymethyl)-3-[2-oxo-2-(3-pyridinylmethylamino)ethyl]-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-6-yl]-4-oxanecarboxamide '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)C(CC(=O)[O-])C=4C=CC=CC4)=O)=O)O)C)(=O)[O-])(=O)[O-])[C@H]([C@H]1O)OP([O-])([O-])=O '
               'NAME: 2-phenylsuccinyl-CoA(5-) REASON: WRONGLY CLASSIFIED '
               'Contains long hydrocarbon chain and multiple amide bonds and '
               'chiral centers indicative of lipopeptide structure\n'
               ' * SMILES: '
               'C[C@@H]1CN(C(=O)C2=C(N=CC(=C2)C3=CN=CN=C3)O[C@H]1CN(C)C(=O)CC4=CC=CC=C4)[C@@H](C)CO '
               'NAME: '
               'N-[[(2R,3R)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-8-(5-pyrimidinyl)-3,4-dihydro-2H-pyrido[2,3-b][1,5]oxazocin-2-yl]methyl]-N-methyl-2-phenylacetamide '
               'REASON: WRONGLY CLASSIFIED Contains long hydrocarbon chain and '
               'multiple amide bonds and chiral centers indicative of '
               'lipopeptide structure\n'
               ' * SMILES: '
               'O([C@H]1[C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H]1OC[C@H]5O[C@@H](O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O[C@H]7[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]7CO)O)[C@@H](O)[C@@H](O[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO[C@]%11(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%11)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]9NC(=O)C)CO)CO)[C@@H]5O)CO)[C@@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@H](O)[C@H](O)[C@H]%13O)CO[C@]%14(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%14)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%12NC(=O)C)CO '
               'NAME: CID 91847795 REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'C[C@@H](N)C(N[C@H](CC=1C=C2C(=CC1)C=CC=C2)C(N[C@@H](C)C(N[C@@H](CC3=CNC4=C3C=CC=C4)C(N[C@H](CC5=CC=CC=C5)C(N[C@@H](CCCCN)C(N)=O)=O)=O)=O)=O)=O '
               'NAME: pralmorelin REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               ' * SMILES: '
               'S(OC[C@H]1OC(OC[C@H]2O[C@H](O)[C@H](NC(=O)C)[C@@H](OC3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO)[C@H]2O)[C@H](NC(=O)C)[C@@H](O)[C@@H]1OC4O[C@@H]([C@H](O)[C@H](O)[C@H]4OC5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)CO)(O)(=O)=O '
               'NAME: '
               '[(2R,3S,4R,5R)-5-Acetamido-6-[[(2R,3R,4R,5R,6S)-5-acetamido-3,6-dihydroxy-4-[(3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]methoxy]-3-[(3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4-hydroxyoxan-2-yl]methyl '
               'hydrogen sulfate REASON: WRONGLY CLASSIFIED Contains long '
               'hydrocarbon chain and multiple amide bonds and chiral centers '
               'indicative of lipopeptide structure\n'
               'False negatives: SMILES: '
               'CCCCCCCC(=O)N(C)[C@H](Cc1ccccc1)C(=O)N[C@@H]1CCCCNC1=O NAME: '
               'ciliatamide B REASON: MISSED Found 2 chiral centers, '
               'indicating insufficient peptide component\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Found 0 amide bonds, need at '
                                               'least 3 for peptide linkage'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Found 0 amide bonds, need at '
                                               'least 3 for peptide linkage'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Found 0 amide bonds, need at '
                                               'least 3 for peptide linkage'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Found 0 amide bonds, need at '
                                               'least 3 for peptide linkage'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Found 0 amide bonds, need at '
                                               'least 3 for peptide linkage'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Found 0 amide bonds, need at '
                                               'least 3 for peptide linkage'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Found 0 amide bonds, need at '
                                               'least 3 for peptide linkage'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'No long hydrocarbon chain '
                                               'found (indicating lipid '
                                               'motif)'},
                                 {   'smiles': '[H][C@@]1(CNc2nc(N)[nH]c(=O)c2N1)[C@@H](O)[C@H](C)O',
                                     'name': 'sapropterin',
                                     'reason': 'Found 0 amide bonds, need at '
                                               'least 3 for peptide linkage'},
                                 {   'smiles': 'CCCCN(CCCC)CCNC(=O)CNC(=O)C1=CC=CS1',
                                     'name': 'N-[2-[2-(dibutylamino)ethylamino]-2-oxoethyl]-2-thiophenecarboxamide',
                                     'reason': 'Found 2 amide bonds, need at '
                                               'least 3 for peptide linkage'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'The analysis of the previous program and its outcomes '
                 'reveals a few areas where potential improvements can be made '
                 'to increase the accuracy of identifying lipopeptide '
                 "structures from their SMILES strings. Let's go through the "
                 'key points:\n'
                 '\n'
                 '1. **Long Hydrocarbon Chain Detection**: The previous '
                 'attempt used a fixed path length of 16 for detecting long '
                 'hydrocarbon chains. However, the definition of "long" can '
                 'vary, especially considering different lipopeptides, which '
                 'may incorporate varying chain lengths. We need to refine '
                 'this by checking for longer or any specific range of path '
                 'lengths that are more representative of lipid chains.\n'
                 '\n'
                 '2. **Amide Bond Detection**: The SMILES patterns used to '
                 'identify amide bonds should capture various peptide linkage '
                 'scenarios. Perhaps, the minimal requirement of two amide '
                 'bonds is too permissive—implementing a condition that checks '
                 'for a higher number might reduce false positives.\n'
                 '\n'
                 '3. **Chiral Centers**: While having several chiral centers '
                 'can imply peptide complexity, it cannot be the sole marker '
                 'for identification since some lipopeptides may have '
                 'structural variations that contain fewer chiral centers. \n'
                 '\n'
                 '4. **Improvement Hypothesis**: Including additional '
                 'structural motifs indicative of lipopeptides such as '
                 'N-terminus (or any other distinctive peptide-related '
                 'functional groups) should improve the specificity of the '
                 'identification process.\n'
                 '\n'
                 "Based on these insights, let's revise the function to "
                 'improve accuracy:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 58,
    'num_false_positives': 100,
    'num_true_negatives': 752,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3670886075949367,
    'recall': 1.0,
    'f1': 0.5370370370370371,
    'accuracy': 0.8901098901098901,
    'negative_predictive_value': 1.0}