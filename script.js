// --- 핵심 계산 로직 (수정 없음) ---
function ellipke(m){let a=1,b=Math.sqrt(1-m),c=1-m/2,d=0;for(;-27>Math.log10(Math.abs(a-b))&&50>d;) {let e=(a+b)/2,f=Math.sqrt(a*b),g=(a-b)/2;c-=Math.pow(2,d)*g*g;a=e;b=f;d++}let e=Math.PI/(2*a);return[e,e*c]}function f24Array(p,a,b,c,d,e,f,g){const h=new Array(p.length),i=a,j=b,k=c,l=d,m=e*e,n=f*f,o=g*g,q=m+o,r=Math.sqrt(q),s=q+n,t=Math.sqrt(s),u=q*s,v=r*t;for(let w=0;w<p.length;w++){const x=p[w],y=Math.sin(x),z=Math.cos(x);let A,B,C,D,E,F;0===r?(B=0,C=-k*Math.sign(f),D=0,E=-j*Math.sign(f),F=l,A=Math.sqrt(j*j+k*k+i*i*z*z-2*i*j*Math.sign(f)*z)):(B=k*g/r,C=-(j*q+k*e*f)/(r*t),D=i*g/t,E=(k*q-j*e*f-l*f*g)/t,F=(l*e-j*g)/r,A=Math.sqrt(j*j+k*k+i*i*((1-n*o/u)*z*z+o/q*y*y+e*f*g/(q*t)*Math.sin(2*x))-(2*i/(v))*(j*e*f-k*q)*z-(2*i*j*g/r)*y));const G=(1+j*j+k*k+i*i+l*l)+2*i*(E*z+F*y),H=4*A/(G+2*A),I=Math.sqrt(H),[,J]=ellipke(H);h[w]=(B*z+C*y+D)*((2-H)*ellipke(H)[0]-2*J)/(I*Math.pow(A,1.5))}return h}function Babic_24(p,a,b,c,d=1e-13){const e=b[0]/p,f=b[1]/p,g=b[2]/p,h=c[0],i=c[1],j=c[2],k=a/p,l=e,m=f,n=g,o=Math.abs(Math.floor(Math.log10(d))),q=Math.pow(2,o-1)+1,r=new Array(q);for(let s=0;s<q;s++)r[s]=2*Math.PI*s/(q-1);const t=f24Array(r,k,l,m,n,h,i,j);let u=new Array(o).fill(0),v=new Array(o).fill(0),w=2*Math.PI;u[0]=w*(t[0]+t[q-1])/2;for(let s=2;s<=o;s++){const x=Math.pow(2,o-s+1);let y=0;for(let z=x/2;z<q-1;z+=x)y+=t[z];v[0]=(u[0]+w*y)/2;for(let z=1;z<s;z++){const A=Math.pow(4,z);v[z]=(A*v[z-1]-u[z-1])/(A-1)}for(let z=0;z<s;z++)u[z]=v[z];w/=2}const s=u[o-1];return 4e-7*Math.PI*a*s}

// --- 전역 변수 ---
let scene, camera, renderer, controls, coilGroup, sweepChartInstance;
let txMeshes = [], rxMeshes = [];
let sweepLabels = [], sweepValues = [];

// --- 테마 및 UI 초기화 ---
document.addEventListener('DOMContentLoaded', () => {
    initTheme();
    initThree();
    initChart();
    
    document.getElementById('calculateBtn').addEventListener('click', calculateAndDisplay);
    document.getElementById('plotBtn').addEventListener('click', plotSweepGraph);
    document.getElementById('saveDataBtn').addEventListener('click', saveSweepData);

    calculateAndDisplay(); // 초기 로드 시 계산 실행
});

function initTheme() {
    const themeToggleBtn = document.getElementById('theme-toggle-btn');
    const sunIcon = `<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/></svg>`;
    const moonIcon = `<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>`;

    const applyTheme = (theme) => {
        document.documentElement.setAttribute('data-theme', theme);
        themeToggleBtn.innerHTML = theme === 'light' ? moonIcon : sunIcon;
        localStorage.setItem('theme', theme);
        updateVisualsTheme(theme);
    };

    themeToggleBtn.addEventListener('click', () => {
        const newTheme = document.documentElement.getAttribute('data-theme') === 'light' ? 'dark' : 'light';
        applyTheme(newTheme);
    });
    
    const savedTheme = localStorage.getItem('theme') || (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    applyTheme(savedTheme);
}

function updateVisualsTheme(theme) {
    const styles = getComputedStyle(document.documentElement);
    const surfaceColor = styles.getPropertyValue('--surface').trim();
    const textColor = styles.getPropertyValue('--text-primary').trim();
    const accentColor = styles.getPropertyValue('--accent').trim();
    const borderColor = styles.getPropertyValue('--border').trim();

    // 3D 뷰 업데이트
    if (scene) {
        scene.background = new THREE.Color(surfaceColor);
    }

    // 차트 업데이트
    if (sweepChartInstance) {
        sweepChartInstance.options.scales.x.ticks.color = textColor;
        sweepChartInstance.options.scales.y.ticks.color = textColor;
        sweepChartInstance.options.scales.x.title.color = textColor;
        sweepChartInstance.options.scales.y.title.color = textColor;
        sweepChartInstance.options.plugins.legend.labels.color = textColor;
        sweepChartInstance.options.scales.x.grid.color = borderColor;
        sweepChartInstance.options.scales.y.grid.color = borderColor;
        sweepChartInstance.data.datasets[0].borderColor = accentColor;
        sweepChartInstance.data.datasets[0].backgroundColor = accentColor + '33'; // Add alpha
        sweepChartInstance.update();
    }
}

// --- 3D 시각화 (Three.js) ---
function initThree() {
    const container = document.getElementById('coilVis');
    scene = new THREE.Scene();
    scene.up.set(0, 0, 1);

    camera = new THREE.PerspectiveCamera(50, container.clientWidth / container.clientHeight, 0.1, 1000);
    camera.up.set(0, 0, 1);
    camera.position.set(20, -20, 20);
    camera.lookAt(0, 0, 0);

    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(container.clientWidth, container.clientHeight);
    container.appendChild(renderer.domElement);

    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;

    const ambientLight = new THREE.AmbientLight(0xffffff, 0.7);
    scene.add(ambientLight);
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
    if (container.clientWidth > 0 && container.clientHeight > 0) {
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
    // 이전 코일 제거
    while(coilGroup.children.length > 0){ 
        coilGroup.remove(coilGroup.children[0]); 
    }

    const inputs = getInputs();
    const tubeRadius = 0.05;

    // 코일 생성 함수
    const createCoil = (params, isTx) => {
        const coil = new THREE.Group();
        const { r_max_cm, N_spiral, p_spiral_cm, N_helix, p_helix_cm } = params;
        
        for (let i = 0; i < N_spiral; i++) {
            const pitch_s = (i < p_spiral_cm.length ? p_spiral_cm[i] : p_spiral_cm.slice(-1)[0]);
            let radius = r_max_cm - p_spiral_cm.slice(0, i).reduce((a, b) => a + b, 0);
            if (radius <= 0) continue;

            for (let j = 0; j < N_helix; j++) {
                const pitch_h = (j < p_helix_cm.length ? p_helix_cm[j] : p_helix_cm.slice(-1)[0]);
                let zPos = p_helix_cm.slice(0, j).reduce((a, b) => a + b, 0);

                const geom = new THREE.TorusGeometry(radius, tubeRadius, 16, 100);
                const mat = new THREE.MeshStandardMaterial({
                    color: isTx ? 0x2563eb : 0xef4444,
                    metalness: 0.5,
                    roughness: 0.5
                });
                const mesh = new THREE.Mesh(geom, mat);
                mesh.position.z = zPos;
                coil.add(mesh);
            }
        }
        return coil;
    };

    // 송신 코일
    const txCoil = createCoil({
        r_max_cm: inputs.r_max_tx_cm, N_spiral: inputs.N_spiral_tx, p_spiral_cm: inputs.p_spiral_tx_cm,
        N_helix: inputs.N_helix_tx, p_helix_cm: inputs.p_helix_tx_cm
    }, true);
    coilGroup.add(txCoil);

    // 수신 코일
    const rxCoil = createCoil({
        r_max_cm: inputs.r_max_rx_cm, N_spiral: inputs.N_spiral_rx, p_spiral_cm: inputs.p_spiral_rx_cm,
        N_helix: inputs.N_helix_rx, p_helix_cm: inputs.p_helix_rx_cm
    }, false);
    
    rxCoil.position.set(inputs.x_off_cm, inputs.y_off_cm, 0);
    rxCoil.rotation.x = inputs.theta;
    rxCoil.rotation.z = inputs.phi; // Tilt around Z-axis first, then X
    coilGroup.add(rxCoil);
}

// --- 차트 (Chart.js) ---
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

// --- 계산 및 플로팅 ---
function getInputs() {
    const p_spiral_tx_cm = document.getElementById('pSpiralTx').value.split(',').map(x => parseFloat(x.trim()) * 1e-1);
    const p_helix_tx_cm = document.getElementById('pHelixTx').value.split(',').map(x => parseFloat(x.trim()) * 1e-1);
    const p_spiral_rx_cm = document.getElementById('pSpiralRx').value.split(',').map(x => parseFloat(x.trim()) * 1e-1);
    const p_helix_rx_cm = document.getElementById('pHelixRx').value.split(',').map(x => parseFloat(x.trim()) * 1e-1);

    return {
        r_max_tx_cm: parseFloat(document.getElementById('rmaxTx').value),
        N_spiral_tx: parseInt(document.getElementById('nSpiralTx').value),
        p_spiral_tx_cm, p_helix_tx_cm,
        N_helix_tx: parseInt(document.getElementById('nHelixTx').value),
        
        r_max_rx_cm: parseFloat(document.getElementById('rmaxRx').value),
        N_spiral_rx: parseInt(document.getElementById('nSpiralRx').value),
        p_spiral_rx_cm, p_helix_rx_cm,
        N_helix_rx: parseInt(document.getElementById('nHelixRx').value),
        
        x_off_cm: parseFloat(document.getElementById('xOffset').value),
        y_off_cm: parseFloat(document.getElementById('yOffset').value),
        z_off_cm: parseFloat(document.getElementById('zOffset').value),
        phi: parseFloat(document.getElementById('phi').value) * Math.PI / 180,
        theta: parseFloat(document.getElementById('theta').value) * Math.PI / 180,
    };
}

function calculateMutualInductance(params) {
    const { r_max_tx_cm, N_spiral_tx, p_spiral_tx_cm, N_helix_tx, p_helix_tx_cm,
            r_max_rx_cm, N_spiral_rx, p_spiral_rx_cm, N_helix_rx, p_helix_rx_cm,
            x_off_cm, y_off_cm, z_off_cm, phi, theta } = params;

    const a = Math.sin(theta) * Math.cos(phi);
    const b = Math.sin(theta) * Math.sin(phi);
    const c = Math.cos(theta);
    const n = [a, b, c];

    const tx_loops = [];
    for (let i = 0; i < N_spiral_tx; i++) {
        const radius = (r_max_tx_cm - p_spiral_tx_cm.slice(0, i).reduce((acc, val) => acc + val, 0)) * 1e-2;
        for (let j = 0; j < N_helix_tx; j++) {
            const zPos = p_helix_tx_cm.slice(0, j).reduce((acc, val) => acc + val, 0) * 1e-2;
            if(radius > 0) tx_loops.push({ r: radius, z: zPos });
        }
    }

    const rx_loops = [];
    for (let i = 0; i < N_spiral_rx; i++) {
        const radius = (r_max_rx_cm - p_spiral_rx_cm.slice(0, i).reduce((acc, val) => acc + val, 0)) * 1e-2;
        for (let j = 0; j < N_helix_rx; j++) {
            const zPos = p_helix_rx_cm.slice(0, j).reduce((acc, val) => acc + val, 0) * 1e-2;
             if(radius > 0) rx_loops.push({ r: radius, z: zPos });
        }
    }

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
    const params = getInputs();
    const M_uH = calculateMutualInductance(params);
    
    document.getElementById('result').innerHTML = `
        <div class="result-item" style="padding: 1rem;">
            <div class="result-label">Mutual Inductance (M)</div>
            <div class="result-value" style="font-size: 1.5rem;">${M_uH.toFixed(6)} μH</div>
        </div>`;

    updateCoilVisuals();
}

async function plotSweepGraph() {
    const progressDiv = document.getElementById('progress');
    const plotBtn = document.getElementById('plotBtn');
    plotBtn.disabled = true;
    plotBtn.classList.add('calculating');

    const sweepVar = document.getElementById('sweepVar').value;
    const sweepMin = parseFloat(document.getElementById('sweepMin').value);
    const sweepMax = parseFloat(document.getElementById('sweepMax').value);
    const sweepStep = parseFloat(document.getElementById('sweepStep').value);

    const points = [];
    for (let v = sweepMin; v <= sweepMax + 1e-9; v += sweepStep) {
        points.push(v);
    }

    sweepLabels = [];
    sweepValues = [];
    
    for (let k = 0; k < points.length; k++) {
        const percent = Math.round(((k + 1) / points.length) * 100);
        progressDiv.textContent = `Calculating... ${percent}%`;
        
        let params = getInputs();
        let val = points[k];

        switch (sweepVar) {
            case 'x': params.x_off_cm = val; break;
            case 'y': params.y_off_cm = val; break;
            case 'z': params.z_off_cm = val; break;
            case 'phi': params.phi = val * Math.PI / 180; break;
            case 'theta': params.theta = val * Math.PI / 180; break;
        }

        const Mval = calculateMutualInductance(params);
        sweepLabels.push(val.toFixed(2));
        sweepValues.push(Mval);

        if (k % 5 === 0 || k === points.length - 1) { // Update chart periodically
            sweepChartInstance.data.labels = sweepLabels;
            sweepChartInstance.data.datasets[0].data = sweepValues;
            sweepChartInstance.options.scales.x.title.text = `${sweepVar} (${sweepVar.length > 2 ? 'deg' : 'cm'})`;
            sweepChartInstance.update('none');
        }
        await new Promise(resolve => setTimeout(resolve, 0)); // Allow UI to update
    }

    progressDiv.textContent = 'Plotting complete.';
    plotBtn.disabled = false;
    plotBtn.classList.remove('calculating');
}

function saveSweepData() {
    if (!sweepLabels.length) {
        alert('No data to save. Please plot first.');
        return;
    }
    const header = `${document.getElementById('sweepVar').value}, M_uH\n`;
    const csvContent = header + sweepLabels.map((label, i) => `${label},${sweepValues[i].toFixed(8)}`).join('\n');
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = "sweep_data.csv";
    link.click();
    URL.revokeObjectURL(link.href);
}
