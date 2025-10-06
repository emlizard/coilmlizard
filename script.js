// --- 1. 핵심 계산 로직 (오류 수정된 버전) ---

// Elliptic integrals via AGM algorithm (원본 코드로 복원)
function ellipke(m) {
    if (m === 1) return [Infinity, 1];
    let a = 1.0;
    let b = Math.sqrt(1 - m);
    let sumE = 1 - m / 2;
    let p = 1;
    while (Math.abs(a - b) > 1e-15) {
        let an = (a + b) / 2.0;
        let bn = Math.sqrt(a * b);
        let cn = (a - b) / 2.0;
        p *= 2;
        sumE -= p * cn * cn / 2;
        a = an;
        b = bn;
    }
    const K = Math.PI / (2 * a);
    const E = K * sumE;
    return [K, E];
}

// f24 integrand evaluation over array of p (원본 코드와 동일)
function f24Array(pArr, alpha, beta, gamma, delta, a, b, c) {
    const results = new Array(pArr.length);
    const h = alpha, e = beta, g = gamma, d = delta;
    const a2 = a * a, b2 = b * b, c2 = c * c;
    const l2 = a2 + c2, l = Math.sqrt(l2);
    const L2 = l2 + b2, L = Math.sqrt(L2);
    const l2L2 = l2 * L2, lL = l * L;

    for (let idx = 0; idx < pArr.length; idx++) {
        const p = pArr[idx];
        const sp = Math.sin(p), cp = Math.cos(p);
        let V, p1, p2, p3, p4, p5;

        if (l < 1e-9) { // Prevent division by zero
            p1 = 0;
            p2 = -g * Math.sign(b);
            p3 = 0;
            p4 = -e * Math.sign(b);
            p5 = d;
            V = Math.sqrt(e*e + g*g + h*h*cp*cp - 2*h*e*Math.sign(b)*cp);
        } else {
            p1 = (g * c) / l;
            p2 = - (e * l2 + g * a * b) / (l * L);
            p3 = (h * c) / L;
            p4 = (g * l2 - e * a * b - d * b * c) / L;
            p5 = (d * a - e * c) / l;
            const term1 = e*e + g*g;
            const term2 = h*h * ( (1 - (b2 * c2) / (l2L2)) * cp*cp + (c2 / l2) * sp*sp + (a * b * c) / (l2 * L) * Math.sin(2*p) );
            const term3 = - (2 * h / (lL)) * (e * a * b - g * l2) * cp;
            const term4 = - (2 * h * e * c / l) * sp;
            V = Math.sqrt(Math.max(0, term1 + term2 + term3 + term4));
        }

        const A = (1 + e*e + g*g + h*h + d*d) + 2*h*(p4 * cp + p5 * sp);
        const m = 4 * V / (A + 2 * V);
        
        if (m > 1 || m < 0 || !isFinite(m)) {
            results[idx] = 0;
            continue;
        }

        const k = Math.sqrt(m);
        const [Kval, Eval] = ellipke(m);
        const PSI = (2 - m) * Kval - 2 * Eval;
        
        if (!isFinite(k) || k < 1e-9 || !isFinite(V) || V < 1e-9) {
             results[idx] = 0;
        } else {
             results[idx] = (p1 * cp + p2 * sp + p3) * (PSI / (k * Math.pow(V, 1.5)));
        }
    }
    return results;
}

// Babic_24 mutual inductance calculation (원본 코드와 동일)
function Babic_24(Rp, Rs, pc, n, tol = 1e-12) {
    const xc = pc[0] / Rp, yc = pc[1] / Rp, zc = pc[2] / Rp;
    const a = n[0], b = n[1], c = n[2];
    const alpha = Rs / Rp, beta = xc, gamma = yc, delta = zc;

    const decdigs = Math.max(8, Math.abs(Math.floor(Math.log10(tol))));
    const nPts = Math.pow(2, decdigs - 1) + 1;
    const pArr = Array.from({length: nPts}, (_, i) => (2 * Math.PI) * i / (nPts - 1));
    
    const romall = f24Array(pArr, alpha, beta, gamma, delta, a, b, c);

    let romPrev = new Array(decdigs).fill(0);
    let romCurr = new Array(decdigs).fill(0);

    let h = 2 * Math.PI;
    romPrev[0] = h * (romall[0] + romall[nPts - 1]) / 2;

    for (let i = 2; i <= decdigs; i++) {
        const step = Math.pow(2, decdigs - i + 1);
        let sumNew = 0;
        for (let j = step / 2; j < nPts - 1; j += step) {
            sumNew += romall[j];
        }
        romCurr[0] = (romPrev[0] + h * sumNew) / 2;
        for (let k = 1; k < i; k++) {
            const factor = Math.pow(4, k);
            romCurr[k] = (factor * romCurr[k - 1] - romPrev[k - 1]) / (factor - 1);
        }
        romPrev = [...romCurr];
        h /= 2;
    }

    const integralVal = romPrev[decdigs - 1];
    const mu0 = 4 * Math.PI * 1e-7;
    return mu0 * Rs * integralVal;
}


// --- 2. UI 및 시각화 로직 (오류 수정 및 개선) ---

// 전역 변수
let scene, camera, renderer, controls, coilGroup, sweepChartInstance;
let sweepLabels = [], sweepValues = [];

// DOM 로드 시 초기화
document.addEventListener('DOMContentLoaded', () => {
    initTheme();
    initThree();
    initChart();
    
    document.getElementById('calculateBtn').addEventListener('click', calculateAndDisplay);
    document.getElementById('plotBtn').addEventListener('click', plotSweepGraph);
    document.getElementById('saveDataBtn').addEventListener('click', saveSweepData);

    calculateAndDisplay(); // 초기 로드 시 계산 및 3D 뷰 생성
});

// 테마 초기화 및 동기화 로직
function initTheme() {
    const themeToggleBtn = document.getElementById('theme-toggle-btn');
    const sunIcon = `<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/></svg>`;
    const moonIcon = `<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>`;

    const applyTheme = (theme) => {
        document.documentElement.setAttribute('data-theme', theme);
        themeToggleBtn.innerHTML = theme === 'light' ? moonIcon : sunIcon;
        localStorage.setItem('theme', theme);
        updateVisualsTheme();
    };

    themeToggleBtn.addEventListener('click', () => {
        const newTheme = document.documentElement.getAttribute('data-theme') === 'light' ? 'dark' : 'light';
        applyTheme(newTheme);
    });
    
    const savedTheme = localStorage.getItem('theme') || (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    applyTheme(savedTheme);
}

function updateVisualsTheme() {
    const styles = getComputedStyle(document.documentElement);
    const surfaceColor = styles.getPropertyValue('--surface').trim();
    const textColor = styles.getPropertyValue('--text-primary').trim();
    const accentColor = styles.getPropertyValue('--accent').trim();
    const borderColor = styles.getPropertyValue('--border').trim();

    if (scene) scene.background = new THREE.Color(surfaceColor);
    if (sweepChartInstance) {
        Object.assign(sweepChartInstance.options.scales.x.ticks, { color: textColor });
        Object.assign(sweepChartInstance.options.scales.y.ticks, { color: textColor });
        Object.assign(sweepChartInstance.options.scales.x.title, { color: textColor });
        Object.assign(sweepChartInstance.options.scales.y.title, { color: textColor });
        Object.assign(sweepChartInstance.options.plugins.legend.labels, { color: textColor });
        Object.assign(sweepChartInstance.options.scales.x.grid, { color: borderColor });
        Object.assign(sweepChartInstance.options.scales.y.grid, { color: borderColor });
        Object.assign(sweepChartInstance.data.datasets[0], { borderColor: accentColor, backgroundColor: accentColor + '33' });
        sweepChartInstance.update();
    }
}

// 3D 시각화 (Three.js)
function initThree() {
    const container = document.getElementById('coilVis');
    scene = new THREE.Scene();
    scene.up.set(0, 0, 1);

    camera = new THREE.PerspectiveCamera(50, container.clientWidth / container.clientHeight, 0.1, 1000);
    camera.position.set(20, -20, 20);
    camera.lookAt(0, 0, 0);

    renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(container.clientWidth, container.clientHeight);
    container.appendChild(renderer.domElement);

    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;

    scene.add(new THREE.AmbientLight(0xffffff, 0.7));
    const dirLight = new THREE.DirectionalLight(0xffffff, 0.5);
    dirLight.position.set(30, -30, 50);
    scene.add(dirLight);

    const axesHelper = new THREE.AxesHelper(10);
    scene.add(axesHelper);

    coilGroup = new THREE.Group();
    scene.add(coilGroup);

    window.addEventListener('resize', onWindowResize);
    animateThree();
}

function onWindowResize() {
    const container = document.getElementById('coilVis');
    if (container.clientWidth > 0) {
        camera.aspect = container.clientWidth / container.clientHeight;
        camera.updateProjectionMatrix();
        renderer.setSize(container.clientWidth, container.clientHeight);
    }
}

function animateThree() {
    requestAnimationFrame(animateThree);
    controls.update();
    renderer.render(scene, camera);
}

function updateCoilVisuals() {
    while(coilGroup.children.length) coilGroup.remove(coilGroup.children[0]);

    const inputs = getInputs();
    const tubeRadius = 0.1;

    const createCoil = (params, isTx) => {
        const { r_max_cm, N_spiral, p_spiral_cm, N_helix, p_helix_cm } = params;
        const group = new THREE.Group();
        
        for (let i = 0; i < N_spiral; i++) {
            const radius = r_max_cm - p_spiral_cm.slice(0, i).reduce((a, b) => a + b, 0);
            if (radius <= 0) continue;
            for (let j = 0; j < N_helix; j++) {
                const zPos = p_helix_cm.slice(0, j).reduce((a, b) => a + b, 0);
                const geom = new THREE.TorusGeometry(radius, tubeRadius, 16, 100);
                const mat = new THREE.MeshStandardMaterial({ color: isTx ? 0x2563eb : 0xef4444, metalness: 0.5, roughness: 0.5 });
                const mesh = new THREE.Mesh(geom, mat);
                mesh.position.z = zPos;
                group.add(mesh);
            }
        }
        return group;
    };

    const txCoil = createCoil({ r_max_cm: inputs.r_max_tx_cm, N_spiral: inputs.N_spiral_tx, p_spiral_cm: inputs.p_spiral_tx_cm, N_helix: inputs.N_helix_tx, p_helix_cm: inputs.p_helix_tx_cm }, true);
    coilGroup.add(txCoil);

    const rxCoil = createCoil({ r_max_cm: inputs.r_max_rx_cm, N_spiral: inputs.N_spiral_rx, p_spiral_cm: inputs.p_spiral_rx_cm, N_helix: inputs.N_helix_rx, p_helix_cm: inputs.p_helix_rx_cm }, false);
    
    // ✅ 3D 위치/방향 로직 수정 (원본 기반)
    const normal = new THREE.Vector3(Math.sin(inputs.theta) * Math.cos(inputs.phi), Math.sin(inputs.theta) * Math.sin(inputs.phi), Math.cos(inputs.theta));
    rxCoil.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), normal);
    rxCoil.position.set(inputs.x_off_cm, inputs.y_off_cm, inputs.z_off_cm);
    coilGroup.add(rxCoil);
}


// 차트 (Chart.js)
function initChart() {
    const ctx = document.getElementById('sweepChart').getContext('2d');
    sweepChartInstance = new Chart(ctx, {
        type: 'line', data: { labels: [], datasets: [{ data: [] }] },
        options: {
            responsive: true, maintainAspectRatio: false,
            scales: { x: { title: { display: true } }, y: { title: { display: true, text: 'Mutual Inductance (μH)' } } },
            plugins: { legend: { display: false } }
        }
    });
}

// 계산 및 플로팅
function getInputs(overrideParams = {}) {
    const getParam = (id, isInt = false, isDeg = false) => {
        const val = parseFloat(document.getElementById(id).value);
        if (isInt) return parseInt(val);
        if (isDeg) return val * Math.PI / 180;
        return val;
    };
    const getPitch = (id) => document.getElementById(id).value.split(',').map(x => parseFloat(x.trim()));

    const params = {
        r_max_tx_cm: getParam('rmaxTx'), N_spiral_tx: getParam('nSpiralTx', true), p_spiral_tx_mm: getPitch('pSpiralTx'),
        N_helix_tx: getParam('nHelixTx', true), p_helix_tx_mm: getPitch('pHelixTx'),
        
        r_max_rx_cm: getParam('rmaxRx'), N_spiral_rx: getParam('nSpiralRx', true), p_spiral_rx_mm: getPitch('pSpiralRx'),
        N_helix_rx: getParam('nHelixRx', true), p_helix_rx_mm: getPitch('pHelixRx'),
        
        x_off_cm: getParam('xOffset'), y_off_cm: getParam('yOffset'), z_off_cm: getParam('zOffset'),
        phi: getParam('phi', false, true), theta: getParam('theta', false, true),
        ...overrideParams
    };
    return params;
}

function calculateMutualInductance(params) {
    const { r_max_tx_cm, N_spiral_tx, p_spiral_tx_mm, N_helix_tx, p_helix_tx_mm, r_max_rx_cm, N_spiral_rx, p_spiral_rx_mm, N_helix_rx, p_helix_rx_mm, x_off_cm, y_off_cm, z_off_cm, phi, theta } = params;

    // ✅ 방향 벡터 계산식 수정 (원본 기반)
    const a = Math.sin(phi) * Math.sin(theta);
    const b = -Math.cos(phi) * Math.sin(theta);
    const c = Math.cos(theta);
    const n = [a, b, c];

    const createLoops = (r_max, Ns, ps, Nh, ph) => {
        const loops = [];
        for (let i = 0; i < Ns; i++) {
            const radius = (r_max - ps.slice(0, i).reduce((acc, val) => acc + val, 0) * 0.1) * 1e-2;
            for (let j = 0; j < Nh; j++) {
                const zPos = (ph.slice(0, j).reduce((acc, val) => acc + val, 0) * 0.1) * 1e-2;
                if(radius > 1e-9) loops.push({ r: radius, z: zPos });
            }
        }
        return loops;
    };

    const tx_loops = createLoops(r_max_tx_cm, N_spiral_tx, p_spiral_tx_mm, N_helix_tx, p_helix_tx_mm);
    const rx_loops = createLoops(r_max_rx_cm, N_spiral_rx, p_spiral_rx_mm, N_helix_rx, p_helix_rx_mm);
    
    let Msum = 0;
    const pc_base = [x_off_cm * 1e-2, y_off_cm * 1e-2, z_off_cm * 1e-2];
    
    for (const tx of tx_loops) {
        for (const rx of rx_loops) {
            const pc = [pc_base[0], pc_base[1], pc_base[2] + rx.z - tx.z];
            let Mij = Babic_24(tx.r, rx.r, pc, n);
            if (!isFinite(Mij)) Mij = 0;
            Msum += Mij;
        }
    }
    return Msum * 1e6; // to μH
}

function calculateAndDisplay() {
    const M_uH = calculateMutualInductance(getInputs());
    document.getElementById('result').innerHTML = `<div class="result-item" style="padding: 1rem;"><div class="result-label">Mutual Inductance (M)</div><div class="result-value" style="font-size: 1.5rem;">${M_uH.toFixed(6)} μH</div></div>`;
    updateCoilVisuals();
}

async function plotSweepGraph() {
    const plotBtn = document.getElementById('plotBtn');
    if (plotBtn.disabled) return;
    plotBtn.disabled = true;
    plotBtn.classList.add('calculating');
    
    const progressDiv = document.getElementById('progress');
    const sweepVar = document.getElementById('sweepVar').value;
    const sweepMin = parseFloat(document.getElementById('sweepMin').value);
    const sweepMax = parseFloat(document.getElementById('sweepMax').value);
    const sweepStep = parseFloat(document.getElementById('sweepStep').value);

    const points = [];
    if (sweepStep > 0) {
        for (let v = sweepMin; v <= sweepMax + 1e-9; v += sweepStep) points.push(v);
    }

    sweepLabels = [];
    sweepValues = [];

    for (let k = 0; k < points.length; k++) {
        const percent = Math.round(((k + 1) / points.length) * 100);
        progressDiv.textContent = `Calculating... ${percent}%`;
        
        let override = {};
        const val = points[k];
        if (sweepVar.length < 3) override[`${sweepVar}_off_cm`] = val;
        else override[sweepVar] = val * Math.PI / 180;

        const Mval = calculateMutualInductance(getInputs(override));
        sweepLabels.push(val.toFixed(2));
        sweepValues.push(Mval);

        if (k % 5 === 0 || k === points.length - 1) {
            sweepChartInstance.data.labels = sweepLabels;
            sweepChartInstance.data.datasets[0].data = sweepValues;
            sweepChartInstance.options.scales.x.title.text = `${sweepVar} (${sweepVar.length > 2 ? 'deg' : 'cm'})`;
            sweepChartInstance.update('none');
        }
        await new Promise(resolve => setTimeout(resolve, 0));
    }

    progressDiv.textContent = 'Plotting complete.';
    plotBtn.disabled = false;
    plotBtn.classList.remove('calculating');
}

function saveSweepData() {
    if (!sweepLabels.length) return alert('No data to save.');
    const header = `${document.getElementById('sweepVar').value}, M_uH\n`;
    const csvContent = header + sweepLabels.map((label, i) => `${label},${sweepValues[i].toFixed(8)}`).join('\n');
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = "sweep_data.csv";
    link.click();
    URL.revokeObjectURL(link.href);
}
