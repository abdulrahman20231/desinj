import math
import pandas as pd
import streamlit as st

def calculate_u(P, T):
    T0 = 459.67
    Tc = 87.89 + T0
    Pc = 1069.87
    a1 = 0.64826
    a2 = 0.12262
    a3 = 0.07287
    a4 = 7.86126
    a5 = 15.68717
    a6 = 1.27533
    a7 = 0.15599
    Mw = 44.1
    A0 = 0.23057
    A1 = 0.00091
    A2 = 0.00221
    A3 = 0.4462
    ur = 0.02179
    y0 = 1.67

    pr = P / Pc
    t = T + T0
    tr = t / Tc

    n = (10 ** (0.3106 - 0.49 * tr + a7 * (tr ** 2))) - 1
    k = (0.62 - 0.23 * tr) + (((a3 / (tr - 0.86)) - 0.037) * pr + ((0.32 * pr ** a4) / (10 ** (a5 * (tr - 1)))))

    po = pr / (tr * (0.1704 * ((tr - 0.7949) ** a1) - 0.0441 * tr - 0.0124))
    rp = (((a2 * po * tr) / pr) - 1) * math.exp(-k * pr)
    rt = po * tr * (0.0162 - 0.0392 * (math.log(tr) ** a6)) * (pr ** n)

    p = po / (1 + rt + rp)
    Z = (P * Mw) / (p * 10.73 * t)

    Eg = (A0 + A1 * t) * p ** y0
    u0 = A2 * t ** A3
    return u0 * math.exp(Eg / t) - ur, Z, p

def get_pseudo(P2, T):
    P1 = 14.7
    pressure_range = list(range(int(P1), int(P2) + 1))
    results = [calculate_u(P, T) for P in pressure_range]
    f_values = [P / (u * Z) for u, Z, P in results]
    pseudo_values = [P1 * f_values[0]]
    
    for i in range(1, len(pressure_range)):
        pseudo = pseudo_values[i - 1] + 2 * ((f_values[i - 1] + f_values[i]) / 2) * (pressure_range[i] - pressure_range[i - 1])
        pseudo_values.append(pseudo)
    
    df = pd.DataFrame({
        'Pressure': pressure_range,
        'u': [result[0] for result in results],
        'Z': [result[1] for result in results],
        'p': [result[2] for result in results],
        'f': f_values,
        'Pseudo': pseudo_values
    })

    return df['Pseudo'].iloc[-1]

def find_pressure_for_pseudo(target_pseudo_value, T, P_min=14.7, P_max=10000, tolerance=1e-6):
    def objective_function(P):
        return get_pseudo(P, T) - target_pseudo_value

    while P_max - P_min > tolerance:
        P_mid = (P_min + P_max) / 2
        if objective_function(P_mid) == 0 or (P_max - P_min) / 2 < tolerance:
            return P_mid
        elif objective_function(P_mid) * objective_function(P_min) < 0:
            P_max = P_mid
        else:
            P_min = P_mid

    return (P_min + P_max) / 2

def calculate_pinj(qinj, pres, T, kh, re, rw):
    pseudo_res = get_pseudo(pres, T)
    ln_term = math.log(re / rw)
    pinj = (pseudo_res) + ((qinj * 1422 * (T + 460) * (ln_term - 0.5)) / (kh))
    pinj = find_pressure_for_pseudo(pinj, T)
    return pinj

# Streamlit app layout and input
st.title('Injection Pressure Calculator for Desired CO2 Injection Rate')
st.header("**Input Parameters:**")
with st.form(key='Parameters'):
    k = st.number_input('Enter Reservoir Permeability (md):', value=300.0, step=0.1)
    h = st.number_input('Enter Reservoir Thickness (ft):', value=300.0, step=0.1)
    kh = k * h
    pres = st.number_input('Enter Reservoir Pressure (Psia):', value=3000.0, step=0.1)
    T = st.number_input('Enter Reservoir Temperature (F):', value=165.0, step=0.1)
    rw = st.number_input('Enter Wellbore Radius (ft):', value=0.333, step=0.001)
    re_max = st.number_input('Enter Reservoir Drainage Radius (ft):', value=1000, step=0.1)
    qinj1_unit = st.selectbox("Select Injection Rate Unit", ["MSCF/D", "ton/day", "Mt/year"], index=0)
    qinj1_value = st.number_input('Enter Injection Rate:', value=3000.0, step=0.1)
    
    # Convert injection rate based on unit selected
    if qinj1_unit == "ton/day":
        qinj1_value *= 34.761636632 / 1000
    elif qinj1_unit == "Mt/year":
        qinj1_value *= 34.761636632 * 1000000 / 365
    
    submitted = st.form_submit_button("Calculate Required Injection Pressure")
    if submitted:
        st.write("Calculating...")  # Debug statement
        pw = calculate_pinj(qinj1_value, pres, T, kh, re_max, rw)
        st.write(f'The Required Injection Pressure (Psia): {pw:.2f}')
