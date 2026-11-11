"""Expand [[CODE]] / [[CODE:text]] link tokens in reference-notebook templates into markdown links.
Usage: python expand.py template.in out.nbsrc   (links are relative to ROOT/reference/)"""
import re, sys
M = {
 "Start": "../00_Start_Here.ipynb",
 "P1.1": "../pass1_big_picture/P1.1_What_MRI_Measures.ipynb",
 "P1.2": "../pass1_big_picture/P1.2_The_MRI_Scanner.ipynb",
 "P1.3": "../pass1_big_picture/P1.3_From_Signal_to_Image.ipynb",
 "P1.4": "../pass1_big_picture/P1.4_Contrast_and_Image_Quality.ipynb",
 "P1.5": "../pass1_big_picture/P1.5_Speed_Artifacts_and_Safety.ipynb",
 "A1": "A1_Spin_Physics", "A2": "A2_Bloch_Equation", "A3": "A3_In_Vivo_Spin_Physics",
 "B1": "B1_MRI_System", "B2": "B2_Magnetic_Fields_in_MRI", "B3": "B3_Main_Magnet_Engineering",
 "B4": "B4_Gradient_System_Engineering", "B5": "B5_RF_System_and_Receive_Chain", "B6": "B6_MRI_Safety",
 "C1": "C1_Pulse_Sequence", "C2": "C2_MRI_Signal_Equation",
 "D1": "D1_Basic_Contrast_Mechanisms", "D2": "D2_In_Vivo_Contrast_Mechanisms",
 "D3": "D3_Gradient_and_Spin_Echoes", "D4": "D4_Image_Contrast_Examples",
 "E1": "E1_RF_Pulses", "E2": "E2_Slice_Selection", "E3": "E3_RF_Pulse_Design",
 "F1": "F1_Spatial_Encoding", "F2": "F2_K-space", "F3": "F3_Image_Reconstruction",
 "G1": "G1_FOV_and_Resolution", "G2": "G2_Signal_to_Noise_Ratio", "G3": "G3_Cartesian_Image_Characteristics",
 "H1": "H1_Volumetric_Imaging", "H2": "H2_Fast_Imaging_Pulse_Sequences", "H3": "H3_Accelerated_Imaging_Methods",
 "I1": "I1_Point_Spread_Function_Analysis", "I2": "I2_Artifacts",
 "J1": "J1_Excitation_K-space", "J2": "J2_Spectral-Spatial_RF_Pulses", "J3": "J3_Fun_with_RF_Pulses",
 "J4": "J4_Advanced_SNR", "J5": "J5_Advanced_Image_Reconstruction",
 "RM1": "RM1_ODEs_for_MRI", "RM2": "RM2_Complex_Numbers_and_Rotations", "RM3": "RM3_Continuous_FT_and_DFT",
 "RM4": "RM4_Fourier_Properties", "RM5": "RM5_Fourier_Pairs_and_Sampling", "RM6": "RM6_2D_Fourier_and_Images",
 "RP1": "RP1_Magnetism_from_Currents", "RP2": "RP2_Induction_and_RF", "RP3": "RP3_Precession_and_Torque",
 "RP4": "RP4_Spin_and_Thermal_Statistics",
 "Notation": "Notation_and_Symbols.ipynb", "Glossary": "Glossary_and_Acronyms.ipynb",
 "Math": "MRI_Math_Concepts.ipynb", "Key": "Key_MRI_Concepts.ipynb", "Questions": "Questions.ipynb",
 "Software": "Software_Resources.ipynb",
}
def path(code):
    p = M[code]
    if p.endswith(".ipynb"):
        return p
    folder = "pass2_deep_dive" if code[0] in "ABCDEFGHIJ" else "refreshers"
    return f"../{folder}/{p}.ipynb"
def sub(m):
    code, _, text = m.group(1).partition(":")
    return f"[{text or code}]({path(code)})"
src = open(sys.argv[1], encoding="utf-8").read()
out = re.sub(r"\[\[([A-Za-z0-9.]+(?::[^\]]*)?)\]\]", sub, src)
left = re.findall(r"\[\[[^\]]*\]\]", out)
if left:
    sys.exit(f"unexpanded tokens: {left[:5]}")
open(sys.argv[2], "w", encoding="utf-8").write(out)
