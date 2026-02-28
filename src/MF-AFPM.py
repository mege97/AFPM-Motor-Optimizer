import femm
import math
import numpy as np
from scipy.optimize import minimize

print("SİSTEM BAŞLATILIYOR: 7 Eksenli High-Fidelity Multiphysics Optimizasyon Motoru Devrede...")


MU_0 = 4 * math.pi * 1e-7
RHO_COPPER_20 = 1.68e-8
H_COERCIVITY = 900000
RPM = 3000
E_MODULUS_STEEL = 200e9
K_H = 0.02; K_E = 1e-4; K_EXC = 1e-5; R_TH = 0.4

DESIGN_PARAMETERS = {
    "phase_current": [1.0, 20.0], "turns_per_coil": [10, 80], "num_poles": [8, 32],
    "outer_diameter": [80.0, 240.0], "inner_diameter": [30.0, 60.0], "num_layers": [2, 4],
    "air_gap": [0.5, 2.0], "magnet_thickness": [2.0, 10.0],
    "rotor_yoke": [3.0, 15.0], "stator_tooth": [10.0, 30.0], "stator_yoke": [3.0, 15.0],       
    "mag_width_ratio": [0.5, 0.85]    
}
def ebob(a, b):
    while b:
        a, b = b, a % b
    return a

def ekok(a, b):
    if a == 0 or b == 0: return 0
    return (a * b) // ebob(a, b)

def mod_1_sacaklanma_kacak(g, h_m, w_slot, w_tooth, w_mag, pole_pitch):
    slot_gap_ratio = w_slot / g
    gamma = (4 / math.pi) * ((slot_gap_ratio / 2) * math.atan(slot_gap_ratio / 2) - math.log(math.sqrt(1 + (slot_gap_ratio / 2)**2)))
    tau_s = w_slot + w_tooth
    k_carter = tau_s / (tau_s - (gamma * g))
    g_eff = g * k_carter
    magnet_gap = pole_pitch - w_mag
    p_gap = w_mag / g_eff
    p_leakage = (h_m / math.pi) * math.log(1 + (math.pi * g_eff / magnet_gap)) if magnet_gap > 0 else 0.0
    k_leakage = p_gap / (p_gap + p_leakage)
    return g_eff, k_leakage

def mod_2_ac_bakir(f_e, d_wire, m_layers, temp_c):
    rho_t = RHO_COPPER_20 * (1 + 0.00393 * (temp_c - 20.0))
    if f_e <= 0: return rho_t, 1.0
    skin_depth = math.sqrt(rho_t / (math.pi * f_e * MU_0))
    h_eq = d_wire * (math.sqrt(math.pi) / 2)
    xi = h_eq / skin_depth
    if xi > 10:
        phi_xi, psi_xi = xi, 2 * xi
    else:
        sinh_2xi, sin_2xi = math.sinh(2 * xi), math.sin(2 * xi)
        cosh_2xi, cos_2xi = math.cosh(2 * xi), math.cos(2 * xi)
        sinh_xi, sin_xi = math.sinh(xi), math.sin(xi)
        cosh_xi, cos_xi = math.cosh(xi), math.cos(xi)
        phi_xi = xi * ((sinh_2xi + sin_2xi) / (cosh_2xi - cos_2xi))
        psi_xi = 2 * xi * ((sinh_xi - sin_xi) / (cosh_xi + cos_xi))
    k_ac = max(1.0, phi_xi + (((m_layers**2) - 1) / 3) * psi_xi)
    return rho_t, k_ac

def mod_3_manyetik_doyum(phi_pole, a_tooth, a_yoke, l_tooth, l_yoke, g_eff, a_gap):
    b_gap = phi_pole / a_gap
    b_tooth = phi_pole / a_tooth
    b_yoke = phi_pole / (2 * a_yoke)
    h_tooth = 14.6 * b_tooth + 31.4 * (b_tooth ** 9.0)
    h_yoke = 14.6 * b_yoke + 31.4 * (b_yoke ** 9.0)
    mmf_gap = (b_gap / MU_0) * g_eff
    mmf_steel = (h_tooth * l_tooth) + (h_yoke * l_yoke)
    k_sat = 1.0 + (mmf_steel / mmf_gap) if mmf_gap > 0 else 1.0
    return b_gap, b_tooth, k_sat

def mod_4_ruzgar_surtunme(rpm, r_out, m_rotor):
    omega = rpm * (2 * math.pi / 60)
    reynolds = (1.225 * omega * (r_out**2)) / 1.81e-5
    c_m = 3.87 / math.sqrt(reynolds) if reynolds < 3e5 else 0.0514 / (reynolds ** 0.2)
    p_windage = 0.5 * c_m * 1.225 * (omega**3) * (r_out**5)
    f_load = m_rotor * 9.81
    p_friction = (0.5 * 0.0015 * f_load * 0.02) * omega
    return p_windage + p_friction

def mod_5_girdap(rpm, num_slots, b_gap, w_slot, g_eff, w_mag, l_mag, h_m):
    f_slot = (rpm / 60.0) * num_slots
    if f_slot <= 0: return 0.0
    slot_gap_ratio = w_slot / g_eff
    beta = 0.5 * (1 - (1 / math.sqrt(1 + (slot_gap_ratio / 2)**2)))
    b_ripple = b_gap * beta  
    v_mag = w_mag * l_mag * h_m
    omega_slot = 2 * math.pi * f_slot
    return (1.0 / 12.0) * 0.625e6 * (omega_slot**2) * (b_ripple**2) * (w_mag**2) * v_mag

def mod_6_esneme(b_gap, r_out, r_in, r_yoke):
    a_active = math.pi * (r_out**2 - r_in**2)
    f_axial = (b_gap**2 * a_active) / (2 * MU_0)
    d_rigidity = (E_MODULUS_STEEL * (r_yoke**3)) / (12 * (1 - 0.3**2))
    q_pressure = f_axial / a_active
    c_factor = (r_out**4 / 64) - ((r_out**2 * r_in**2) / 8) + (r_in**4 / 64) + ((r_out**2 * r_in**2) / 16)*math.log(r_out/r_in)
    deflection = (q_pressure / d_rigidity) * c_factor if d_rigidity > 0 else float('inf')
    return deflection

def mod_7_vuruntu(num_poles, num_slots, b_gap, r_out, r_in, w_slot, g_eff):
    l_c_m = ekok(int(num_poles), int(num_slots))
    if l_c_m == 0: return 0.0
    a_active = math.pi * (r_out**2 - r_in**2)
    cogging_factor = (math.pi / 4) * (w_slot / g_eff) * math.exp(-0.1 * (l_c_m / num_poles))
    return ((b_gap**2 * a_active * g_eff) / (2 * MU_0)) * cogging_factor

def calculate_physics(design):

    r_out = design["outer_diameter"] / 2000; r_in = design["inner_diameter"] / 2000
    r_avg = (r_out + r_in) / 2; l_mag = r_out - r_in
    g = design["air_gap"] / 1000; h_m = design["magnet_thickness"] / 1000
    r_yoke = design["rotor_yoke"] / 1000; s_tooth = design["stator_tooth"] / 1000; s_yoke = design["stator_yoke"] / 1000
    p = int(design["num_poles"]); num_slots = int(p * 1.5)
    
    pole_pitch = (2 * math.pi * r_avg) / p
    slot_pitch = (2 * math.pi * r_avg) / num_slots
    w_tooth = slot_pitch * 0.5; w_slot = slot_pitch * 0.5
    w_mag = pole_pitch * design["mag_width_ratio"]
    
    g_eff, k_leakage = mod_1_sacaklanma_kacak(g, h_m, w_slot, w_tooth, w_mag, pole_pitch)
    
    a_gap = w_mag * l_mag
    phi_ideal = (H_COERCIVITY * h_m) / (g_eff / (MU_0 * a_gap))
    phi_real = phi_ideal * k_leakage
    b_gap, b_tooth, k_sat = mod_3_manyetik_doyum(phi_real, w_tooth*l_mag, s_yoke*l_mag, s_tooth, pole_pitch, g_eff, a_gap)
    f_e = (p / 2) * (RPM / 60)
    wire_area = (w_slot * s_tooth * 0.4) / design["turns_per_coil"]
    d_wire = math.sqrt((4 * wire_area) / math.pi)
    rho_t, k_ac = mod_2_ac_bakir(f_e, d_wire, design["num_layers"], temp_c=80.0)
    
    p_m_coil = 2 * l_mag + 2 * slot_pitch
    r_dc = rho_t * (p_m_coil * design["turns_per_coil"]) / wire_area
    r_ac = r_dc * k_ac
    p_cu_ac = 3 * (design["phase_current"]**2) * r_ac * num_slots
    
    m_rotor = math.pi * (r_out**2) * r_yoke * 7850 # Çelik yoğunluğu
    p_mech = mod_4_ruzgar_surtunme(RPM, r_out, m_rotor)
    
    p_eddy = mod_5_girdap(RPM, num_slots, b_gap, w_slot, g_eff, w_mag, l_mag, h_m) * p
    
    vol_iron = math.pi * (r_out**2 - r_in**2) * (r_yoke + s_tooth + s_yoke)
    p_core = ((K_H * f_e * b_tooth**2) + (K_E * f_e**2 * b_tooth**2)) * vol_iron * 7850
    
    tq_em = 1.5 * (p / 2) * phi_real * design["phase_current"]
    
    t_cogging = mod_7_vuruntu(p, num_slots, b_gap, r_out, r_in, w_slot, g_eff)
    
    deflection = mod_6_esneme(b_gap, r_out, r_in, r_yoke)
    
    p_out = (tq_em * RPM * 2 * math.pi / 60) - p_mech
    p_total_loss = p_cu_ac + p_core + p_mech + p_eddy
    efficiency = p_out / (p_out + p_total_loss + 1e-6)
    
    return {
        "torque": tq_em, "efficiency": efficiency, "cogging": t_cogging, 
        "deflection": deflection, "g_limit": g/2, "k_sat": k_sat,
        "temp": 40 + (p_total_loss * R_TH)
    }
def optimize():
    keys = list(DESIGN_PARAMETERS.keys()); bounds = [DESIGN_PARAMETERS[k] for k in keys]; x0 = [sum(b)/2 for b in bounds]
    cons = [
        {'type': 'ineq', 'fun': lambda x: calculate_physics({keys[i]: x[i] for i in range(len(keys))})["efficiency"] - 0.90}, # %90 Verim zorunlu
        {'type': 'ineq', 'fun': lambda x: 120.0 - calculate_physics({keys[i]: x[i] for i in range(len(keys))})["temp"]},      # 120 derece termal limit
        {'type': 'ineq', 'fun': lambda x: calculate_physics({keys[i]: x[i] for i in range(len(keys))})["g_limit"] - calculate_physics({keys[i]: x[i] for i in range(len(keys))})["deflection"]}, # Çarpışma engelleme
        {'type': 'ineq', 'fun': lambda x: 1.5 - calculate_physics({keys[i]: x[i] for i in range(len(keys))})["k_sat"]}        # Doyum limiti (1.5 max)
    ]
    res = minimize(lambda x, k: -calculate_physics({k[i]: x[i] for i in range(len(k))})["torque"], x0, args=(keys,), method='SLSQP', bounds=bounds, constraints=cons)
    return {keys[i]: res.x[i] for i in range(len(keys))}

print("Milyonlarca kombinasyon test ediliyor... (Bu işlem ağır fizik kuralları içerdiğinden 10-20 saniye sürebilir)")
opt = optimize()
res = calculate_physics(opt)

print(f"\n[BÜYÜK BİRLEŞME TAMAMLANDI - OPTİMUM DEĞERLER]")
print(f"Net Üretilen Tork: {res['torque']:.2f} Nm | Verim: %{res['efficiency']*100:.1f}")
print(f"Vuruntu Torku (Cogging): {res['cogging']:.3f} Nm | Rotor Esnemesi: {res['deflection']*1000:.3f} mm")

try:
    print("\nSonuçlar sonlu elemanlar laboratuvarına (FEMM) aktarılıyor...")
    femm.openfemm()
    femm.newdocument(0)
    femm.mi_probdef(0, 'millimeters', 'planar', 1e-8, 0, 30)
    
    femm.mi_getmaterial('Air')
    femm.mi_getmaterial('M-19 Steel')
    femm.mi_getmaterial('NdFeB 52 MGOe')
    femm.mi_getmaterial('Copper')
    
    g = opt["air_gap"]; h_m = opt["magnet_thickness"]
    r_yoke = opt["rotor_yoke"]; s_tooth = opt["stator_tooth"]; s_yoke = opt["stator_yoke"]
    p = int(opt["num_poles"])
    r_avg = (opt["outer_diameter"] + opt["inner_diameter"]) / 4
    pole_pitch = (2 * math.pi * r_avg) / p
    
    y_rotor_alt = 0.0; y_rotor_ust = y_rotor_alt + r_yoke
    y_magnet_alt = y_rotor_ust; y_magnet_ust = y_magnet_alt + h_m  
    y_stator_alt = y_magnet_ust + g; y_stator_ust = y_stator_alt + s_tooth
    y_yoke_ust = y_stator_ust + s_yoke

    def draw_rect(x1, y1, x2, y2):
        femm.mi_drawline(float(x1), float(y1), float(x2), float(y1))
        femm.mi_drawline(float(x2), float(y1), float(x2), float(y2))
        femm.mi_drawline(float(x2), float(y2), float(x1), float(y2))
        femm.mi_drawline(float(x1), float(y2), float(x1), float(y1))

    def assign_mat(x, y, mat_name, mag_dir=0):
        femm.mi_addblocklabel(float(x), float(y))
        femm.mi_selectlabel(float(x), float(y))
        femm.mi_setblockprop(mat_name, 1, 0, '<None>', mag_dir, 0, 0)
        femm.mi_clearselected()

    draw_rect(-5, -5, (pole_pitch * 4) + 5, y_yoke_ust + 5)
    assign_mat(-2, -2, 'Air')
    assign_mat((pole_pitch * 4) / 2, y_magnet_ust + (g / 2), 'Air') 
    
    draw_rect(0, y_rotor_alt, pole_pitch * 4, y_rotor_ust)
    assign_mat((pole_pitch * 4) / 2, y_rotor_alt + (r_yoke / 2), 'M-19 Steel')
    
    w_mag = pole_pitch * opt["mag_width_ratio"]
    for i in range(4):
        x_start = (i * pole_pitch) + ((pole_pitch - w_mag) / 2)
        draw_rect(x_start, y_magnet_alt, x_start + w_mag, y_magnet_ust)
        assign_mat(x_start + (w_mag / 2), y_magnet_alt + (h_m / 2), 'NdFeB 52 MGOe', 90 if i % 2 == 0 else -90)
        
    for i in range(4):
        x_tooth = (i * pole_pitch) + (pole_pitch * 0.25)
        w_tooth = pole_pitch * 0.5
        draw_rect(x_tooth, y_stator_alt, x_tooth + w_tooth, y_stator_ust)
        assign_mat(x_tooth + (w_tooth / 2), y_stator_alt + (s_tooth / 2), 'M-19 Steel')
        
        x_slot = x_tooth + w_tooth
        if i < 3: 
            draw_rect(x_slot, y_stator_alt, x_slot + w_tooth, y_stator_ust)
            assign_mat(x_slot + (w_tooth / 2), y_stator_alt + (s_tooth / 2), 'Copper')

    draw_rect(0, y_stator_ust, pole_pitch * 4, y_yoke_ust)
    assign_mat((pole_pitch * 4) / 2, y_stator_ust + (s_yoke / 2), 'M-19 Steel')
    femm.mi_zoomnatural()
    print("\nAŞAMA 3: Analiz (Maxwell Denklemleri) Çözülüyor...")
    
    import os
    kayit_yolu = os.path.join(os.getcwd(), "Eksenel_Motor_Optimizasyon.fem")
    femm.mi_saveas(kayit_yolu)
    femm.mi_analyze(1)
    femm.mi_loadsolution()
    femm.mo_showdensityplot(1, 0, 2.0, 0, 'bmag')
    femm.mo_zoomnatural()
    
    print("\n[MUAZZAM BAŞARI] Nihai analiz tamamlandı! Lütfen FEMM Post-Processor ekranına geçin.")

except Exception as e:
    print(f"\n[HATA]: {e}")

input("\nProgramı kapatmak için ENTER'a basın...")