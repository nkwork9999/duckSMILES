use ducksmiles_analytics::execute as run;
fn number(j: &str, k: &str) -> f64 {
    let s = j.find(&format!("\"{k}\":")).unwrap() + k.len() + 3;
    j[s..]
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
fn values(j: &str, k: &str) -> Vec<f64> {
    let s = j.find(&format!("\"{k}\":[")).unwrap() + k.len() + 4;
    let t = &j[s..];
    let e = t.find(']').unwrap();
    if e == 0 {
        vec![]
    } else {
        t[..e]
            .split(',')
            .map(|x| x.trim().parse().unwrap())
            .collect()
    }
}
fn matrix(j: &str, k: &str) -> Vec<Vec<f64>> {
    let s = j.find(&format!("\"{k}\":[")).unwrap() + k.len() + 4;
    let t = &j[s..];
    let e = t.find("]]").unwrap() + 1;
    t[..e]
        .split("],")
        .map(|s| {
            s.trim()
                .trim_matches(['[', ']'])
                .split(',')
                .map(|x| x.trim().parse().unwrap())
                .collect()
        })
        .collect()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-8 * (1. + b.abs()), "{a} != {b}")
}
fn check(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        close(*x, *y)
    }
}

#[test]
fn gaussian_rts_matches_joint_normal_conditioning_and_observation_forecast_variance() {
    let j = run(
        "kalman_smoother",
        &[vec![1., 1.], vec![2., 1.], vec![1.5, 1.]],
        &[1., 2., 0.95, 0., 1., 1., 0.2, 1., 0.5],
    )
    .unwrap();
    let m = matrix(&j, "smoothed_means");
    check(
        &m.iter().map(|x| x[0]).collect::<Vec<_>>(),
        &[1.1212121212121213, 1.393939393939394, 1.4242424242424243],
    );
    check(
        &values(&j, "forecast_means"),
        &[1.4242424242424243, 1.4242424242424243],
    );
    check(
        &values(&j, "forecast_observation_variances"),
        &[0.9390572390572391, 1.1390572390572393],
    );
    close(number(&j, "log_likelihood"), -4.159068440286804);
    let mut args = vec![2., 0.95];
    args.extend(values(&j, "prediction_model"));
    let p = run("state_space_predict", &[], &args).unwrap();
    check(&values(&p, "forecast_means"), &values(&j, "forecast_means"));
    check(
        &values(&p, "forecast_observation_variances"),
        &values(&j, "forecast_observation_variances"),
    );
}
#[test]
fn missing_observation_skips_update_but_retains_elapsed_process_variance() {
    let a = [1., 1., 0.95, 0., 1., 1., 0.2, 1., 0.5];
    let j = run(
        "kalman_smoother",
        &[vec![1., 1.], vec![999., 0.], vec![1.5, 1.]],
        &a,
    )
    .unwrap();
    let k = run(
        "kalman_smoother",
        &[vec![1., 1.], vec![-999., 0.], vec![1.5, 1.]],
        &a,
    )
    .unwrap();
    check(&values(&j, "forecast_means"), &values(&k, "forecast_means"));
    close(number(&j, "observed_count"), 2.);
    let p = run("kalman_smoother", &[vec![0., 0.]], &a).unwrap();
    check(&values(&p, "forecast_observation_variances"), &[1.7]);
}
#[test]
fn invalid_covariance_or_incoherent_saved_model_is_rejected() {
    assert!(run(
        "kalman_smoother",
        &[vec![0., 1.]],
        &[1., 1., 0.95, 0., 1., 1., -0.1, 1., 0.5]
    )
    .is_err());
    assert!(run(
        "kalman_smoother",
        &[vec![0., 1.]],
        &[1., 1., 0.95, 0., 1., 1., 0.1, 1., 0.]
    )
    .is_err());
    assert!(run("state_space_predict", &[], &[1., 0.95, 1., 0., 1000000.]).is_err());
}

#[test]
fn additive_ets_six_families_match_independent_dense_innovations_profile() {
    {
        let rows = vec![
            vec![2.0000342874782],
            vec![1.8167048890869033],
            vec![2.166233669768416],
            vec![2.4319014093536797],
            vec![1.9145671885706752],
            vec![2.2319058081630043],
            vec![1.8208108446820492],
            vec![1.9874718189737965],
            vec![2.2978627440436132],
            vec![2.171892731841144],
            vec![1.5485318747658159],
            vec![2.371961926059055],
            vec![1.5798888805582152],
            vec![1.5054699738486075],
            vec![1.8688839526119838],
            vec![1.8693906291991589],
            vec![1.6827784209059515],
            vec![1.6436353611943137],
            vec![1.4872856545940125],
            vec![1.7898123463522735],
            vec![1.7817886754290673],
            vec![1.7298807419968385],
            vec![1.3939250911178402],
            vec![1.9804510339191803],
            vec![1.718139940995836],
            vec![2.00374632849518],
            vec![1.9397600805686943],
            vec![1.9440245159758172],
            vec![1.8919822230680545],
            vec![1.983588592247076],
            vec![1.8695526314740125],
            vec![1.8452449541305507],
            vec![2.27824081564765],
            vec![1.4402425369556728],
            vec![2.089198068409359],
            vec![1.7859874498349018],
            vec![1.84274268915339],
            vec![1.8408527381652366],
            vec![2.4155260626863333],
            vec![1.9704541634733614],
            vec![2.0758589539052767],
            vec![1.9784078013597042],
            vec![1.7903550166297706],
            vec![2.4187116629822465],
            vec![2.5477971654466813],
            vec![2.5818108683295393],
            vec![2.2634909206544567],
            vec![2.864083210376393],
            vec![2.1230315818103582],
            vec![2.503111216236104],
            vec![2.707902295132613],
            vec![2.380642221081157],
            vec![2.6518033667795993],
            vec![2.3657675398169404],
            vec![2.2822055282847358],
            vec![2.6210995515592215],
            vec![2.8963978093107046],
            vec![2.442351394553133],
            vec![2.9629388974197512],
            vec![2.634222261380734],
            vec![2.944554092838275],
            vec![2.2778255027858134],
            vec![2.2662486809267115],
            vec![2.8731271073756504],
        ];
        let j = run(
            "state_space_smoothing",
            &rows,
            &[0., 4., 5., 0.95, 2500., 1e-8],
        )
        .unwrap();
        assert!(
            (number(&j, "sigma_squared") - 0.07810150043837706).abs() < 2e-6,
            "{j}"
        );
        let par = values(&j, "smoothing_parameters");
        let want = [0.25402264281868775, 0.0, 0.0, 1.0];
        for i in 0..4 {
            assert!((par[i] - want[i]).abs() < 0.0004, "{par:?} != {want:?}");
        }
        assert!(j.contains("\"converged\":true"), "{j}");
        let mut args = vec![5., 0.95];
        args.extend(values(&j, "prediction_model"));
        let k = run("state_space_predict", &[], &args).unwrap();
        check(&values(&j, "forecast_means"), &values(&k, "forecast_means"));
        check(
            &values(&j, "forecast_observation_variances"),
            &values(&k, "forecast_observation_variances"),
        );
    }
    {
        let rows = vec![
            vec![1.6626097081966815],
            vec![1.7705840713758443],
            vec![1.9397454584868041],
            vec![1.754161324396973],
            vec![2.0539106436982744],
            vec![1.2819683522277132],
            vec![2.0591761506690203],
            vec![2.0248014971805706],
            vec![2.039623569722697],
            vec![2.021372851496684],
            vec![1.8891996542667686],
            vec![1.747561605910617],
            vec![2.6317519147988158],
            vec![2.0810541556543023],
            vec![1.9346674933455021],
            vec![1.842963369785557],
            vec![2.476767661913438],
            vec![2.4018623473471643],
            vec![2.2025159146465723],
            vec![2.2196515612117476],
            vec![1.9885555278059235],
            vec![2.3814546142367656],
            vec![2.681549374962123],
            vec![2.093418161658853],
            vec![1.9993813908943137],
            vec![2.2931832453062593],
            vec![2.0005311581650376],
            vec![2.1446567883306544],
            vec![2.060942777860613],
            vec![1.9924367961595384],
            vec![1.7426422650927451],
            vec![1.7487811508188398],
            vec![1.8667666896267678],
            vec![1.8749630386571063],
            vec![1.5147894966143718],
            vec![1.7641535006572688],
            vec![2.163783199858959],
            vec![1.9055827416659545],
            vec![1.815000694487582],
            vec![2.042234388344521],
            vec![1.5488171730504674],
            vec![1.3242133854602947],
            vec![1.601660753502369],
            vec![1.742859919842895],
            vec![0.9395278118702701],
            vec![1.2503170101122227],
            vec![0.9030492470867388],
            vec![1.1640701356786591],
            vec![1.0648860813521794],
            vec![1.0272040546680563],
            vec![0.7389537206451986],
            vec![0.5498415983641655],
            vec![0.9916866616349631],
            vec![0.6856048478605895],
            vec![0.41184776759606845],
            vec![0.631453501634242],
            vec![-0.1343434483661825],
            vec![0.2783277130478663],
            vec![-0.24271369475340235],
            vec![-0.43762067390722664],
            vec![-0.8220457190605742],
            vec![-0.6387083572745704],
            vec![-0.6710628814759542],
            vec![-1.2071796774739205],
        ];
        let j = run(
            "state_space_smoothing",
            &rows,
            &[1., 4., 5., 0.95, 2500., 1e-8],
        )
        .unwrap();
        assert!(
            (number(&j, "sigma_squared") - 0.07635640195659274).abs() < 2e-6,
            "{j}"
        );
        let par = values(&j, "smoothing_parameters");
        let want = [0.23592000744415562, 0.054969774613521184, 0.0, 1.0];
        for i in 0..4 {
            assert!((par[i] - want[i]).abs() < 0.0004, "{par:?} != {want:?}");
        }
        assert!(j.contains("\"converged\":true"), "{j}");
        let mut args = vec![5., 0.95];
        args.extend(values(&j, "prediction_model"));
        let k = run("state_space_predict", &[], &args).unwrap();
        check(&values(&j, "forecast_means"), &values(&k, "forecast_means"));
        check(
            &values(&j, "forecast_observation_variances"),
            &values(&k, "forecast_observation_variances"),
        );
    }
    {
        let rows = vec![
            vec![2.325929616024174],
            vec![2.4575581890109217],
            vec![2.4928945155338624],
            vec![1.7172336695987118],
            vec![3.274703855320497],
            vec![3.2659580230466276],
            vec![2.952373500922584],
            vec![3.28204164838051],
            vec![3.3867122953541124],
            vec![3.133030669652646],
            vec![3.859723118130578],
            vec![4.219524734906579],
            vec![4.088024432596456],
            vec![4.587465508412508],
            vec![4.60168707610557],
            vec![5.1882575208230755],
            vec![5.488684014396197],
            vec![5.717088216001393],
            vec![6.106745469461851],
            vec![5.8596174915308215],
            vec![6.400173260302727],
            vec![6.665637263519948],
            vec![6.535110700801465],
            vec![7.229073818320764],
            vec![6.978918934063053],
            vec![7.343670455687638],
            vec![7.7826499748212745],
            vec![7.849890351395448],
            vec![7.592298629241025],
            vec![8.22637055362785],
            vec![8.511356671075276],
            vec![8.420131627045519],
            vec![7.9036228936183335],
            vec![8.702960415931182],
            vec![8.509316704802526],
            vec![8.269743613962437],
            vec![8.0352357209815],
            vec![8.44100835562223],
            vec![8.232395719482453],
            vec![8.146680290755846],
            vec![8.671052543002558],
            vec![8.312535444070388],
            vec![8.480817457265498],
            vec![8.131424669169775],
            vec![8.172624769706706],
            vec![8.009871352755601],
            vec![8.505921912900433],
            vec![8.280648175670008],
            vec![8.214248328300748],
            vec![8.302982859740974],
            vec![8.201694237373985],
            vec![8.226815921630843],
            vec![8.44546801231217],
            vec![7.9511832197109324],
            vec![8.22209366587337],
            vec![8.049613230803793],
            vec![7.668003308987102],
            vec![7.5003942280954465],
            vec![7.581975324565266],
            vec![7.382755953626837],
            vec![7.324288210902362],
            vec![7.07056676922803],
            vec![6.563098537425927],
            vec![6.744555051968641],
        ];
        let j = run(
            "state_space_smoothing",
            &rows,
            &[2., 4., 5., 0.95, 2500., 1e-8],
        )
        .unwrap();
        assert!(
            (number(&j, "sigma_squared") - 0.09517398282650381).abs() < 2e-6,
            "{j}"
        );
        let par = values(&j, "smoothing_parameters");
        let want = [
            0.21646197231195444,
            0.14644177472494185,
            0.0,
            0.9688540005418151,
        ];
        for i in 0..4 {
            assert!((par[i] - want[i]).abs() < 0.0004, "{par:?} != {want:?}");
        }
        assert!(j.contains("\"converged\":true"), "{j}");
        let mut args = vec![5., 0.95];
        args.extend(values(&j, "prediction_model"));
        let k = run("state_space_predict", &[], &args).unwrap();
        check(&values(&j, "forecast_means"), &values(&k, "forecast_means"));
        check(
            &values(&j, "forecast_observation_variances"),
            &values(&k, "forecast_observation_variances"),
        );
    }
    {
        let rows = vec![
            vec![2.541740239185275],
            vec![1.1354716686692643],
            vec![1.8133426125911323],
            vec![2.984553166859041],
            vec![2.2381200961088865],
            vec![1.7667538420900875],
            vec![1.8083769229688151],
            vec![2.2820932904597027],
            vec![2.2392371174021286],
            vec![1.384380885890391],
            vec![1.7617314483106719],
            vec![2.1136919981922437],
            vec![2.4038756167034006],
            vec![1.6650165932745717],
            vec![1.88627959567044],
            vec![2.7856748547465897],
            vec![2.3793618213413406],
            vec![1.5023069087808292],
            vec![2.285133300172262],
            vec![2.358244298003024],
            vec![2.567862595725126],
            vec![1.2910755038764428],
            vec![2.1112132491239146],
            vec![2.352215340533903],
            vec![2.4857994796850607],
            vec![1.2517617203156002],
            vec![1.9116831621233241],
            vec![1.9208235647835024],
            vec![2.1628497213009314],
            vec![1.5453057251888729],
            vec![1.901227588798789],
            vec![2.5411273946099247],
            vec![2.47499349209667],
            vec![1.5806140620875533],
            vec![2.3885037340274327],
            vec![2.0301863571486902],
            vec![2.4433510098600704],
            vec![1.389361172646904],
            vec![2.2360703470618915],
            vec![2.2237853309044295],
            vec![2.023981850096428],
            vec![1.3830972798783199],
            vec![1.9485073701276494],
            vec![1.8626929612103842],
            vec![1.6548174888396692],
            vec![0.7646989692069712],
            vec![1.2508197534200638],
            vec![2.017964498121379],
            vec![1.9010526826388197],
            vec![1.0781082046028576],
            vec![1.9534694131651658],
            vec![2.2967817420998884],
            vec![2.5834327911290123],
            vec![1.7726051254277806],
            vec![2.1961616446842203],
            vec![2.6967812965354923],
            vec![2.592734161680808],
            vec![1.6443334657497544],
            vec![2.406191465203994],
            vec![2.7088576519251673],
            vec![2.270284743765771],
            vec![1.5086305619185816],
            vec![2.3560995731869405],
            vec![2.1739938115098676],
            vec![1.7756944269074149],
            vec![1.501735521387791],
            vec![1.9847612802817167],
            vec![1.733879901328872],
            vec![1.777216359472379],
            vec![1.589075727010331],
            vec![1.80823614869837],
            vec![2.0177932582030427],
            vec![2.0398770646692124],
            vec![1.9944100506738023],
            vec![2.641644129785287],
            vec![2.3112864395683324],
            vec![2.553118322276343],
            vec![1.9948073826716062],
            vec![2.5769023461396268],
            vec![1.8046589911493545],
        ];
        let j = run(
            "state_space_smoothing",
            &rows,
            &[3., 4., 5., 0.95, 2500., 1e-8],
        )
        .unwrap();
        assert!(
            (number(&j, "sigma_squared") - 0.08082745255915653).abs() < 2e-6,
            "{j}"
        );
        let par = values(&j, "smoothing_parameters");
        let want = [0.39302034770016303, 0.0, 0.30652329766249403, 1.0];
        for i in 0..4 {
            assert!((par[i] - want[i]).abs() < 0.0004, "{par:?} != {want:?}");
        }
        assert!(j.contains("\"converged\":true"), "{j}");
        let mut args = vec![5., 0.95];
        args.extend(values(&j, "prediction_model"));
        let k = run("state_space_predict", &[], &args).unwrap();
        check(&values(&j, "forecast_means"), &values(&k, "forecast_means"));
        check(
            &values(&j, "forecast_observation_variances"),
            &values(&k, "forecast_observation_variances"),
        );
    }
    {
        let rows = vec![
            vec![1.881310109692],
            vec![1.425871828200431],
            vec![1.8367028789239597],
            vec![1.947416378919932],
            vec![2.1158252171227625],
            vec![1.4185310848697779],
            vec![1.686404724978314],
            vec![2.5516127598397635],
            vec![1.5293860258745502],
            vec![1.351767507293837],
            vec![1.3766349387979646],
            vec![2.336201711064073],
            vec![1.280830889226631],
            vec![1.1754283230235874],
            vec![1.555503810156635],
            vec![2.187550973822845],
            vec![1.7872997544188232],
            vec![1.395407990634653],
            vec![1.3759494851552958],
            vec![1.492571364767481],
            vec![1.720132477054566],
            vec![1.0999156459085686],
            vec![1.9193852697165044],
            vec![2.315032589717963],
            vec![1.942177467137486],
            vec![1.8689398725547337],
            vec![1.9545477412782282],
            vec![2.7113148009969157],
            vec![2.2082708113648346],
            vec![1.615766493192533],
            vec![2.1598098219405504],
            vec![2.8653457302609],
            vec![2.7804893876565755],
            vec![1.9542530868268586],
            vec![2.6184401671035937],
            vec![3.198144999359686],
            vec![3.2080337666891],
            vec![2.806669794187212],
            vec![2.9517554118915745],
            vec![4.228541146698178],
            vec![4.452922849502518],
            vec![3.6135261068643523],
            vec![4.357758914845486],
            vec![5.518658806910794],
            vec![5.540437439545374],
            vec![4.954850829358626],
            vec![5.565068944418184],
            vec![6.538309836397082],
            vec![6.447721771072628],
            vec![6.41636070544924],
            vec![6.902950019956489],
            vec![7.922055407378174],
            vec![7.791424257475039],
            vec![7.519523951747322],
            vec![8.286631950056748],
            vec![9.091304836169048],
            vec![8.975207804806303],
            vec![9.125380419777986],
            vec![9.481293090490372],
            vec![10.802104883264798],
            vec![10.263570286231436],
            vec![10.162421920544329],
            vec![10.871833985217549],
            vec![12.179014271453166],
            vec![11.723612653276657],
            vec![11.579157600180185],
            vec![12.185920159712278],
            vec![13.00850829564617],
            vec![12.925588126190286],
            vec![12.586201268346922],
            vec![13.768974640930674],
            vec![14.472436161855114],
            vec![14.29001933186793],
            vec![13.741527355427067],
            vec![14.542739890147825],
            vec![16.11156156862649],
            vec![15.445822818880096],
            vec![15.537956540755745],
            vec![15.708507387494933],
            vec![16.88077332818098],
        ];
        let j = run(
            "state_space_smoothing",
            &rows,
            &[4., 4., 5., 0.95, 2500., 1e-8],
        )
        .unwrap();
        assert!(
            (number(&j, "sigma_squared") - 0.06688307649098749).abs() < 2e-6,
            "{j}"
        );
        let par = values(&j, "smoothing_parameters");
        let want = [
            0.2868909371187082,
            0.10990493673439801,
            1.427974375528223e-09,
            1.0,
        ];
        for i in 0..4 {
            assert!((par[i] - want[i]).abs() < 0.0004, "{par:?} != {want:?}");
        }
        assert!(j.contains("\"converged\":true"), "{j}");
        let mut args = vec![5., 0.95];
        args.extend(values(&j, "prediction_model"));
        let k = run("state_space_predict", &[], &args).unwrap();
        check(&values(&j, "forecast_means"), &values(&k, "forecast_means"));
        check(
            &values(&j, "forecast_observation_variances"),
            &values(&k, "forecast_observation_variances"),
        );
    }
    {
        let rows = vec![
            vec![2.604190208411519],
            vec![1.6266917645665346],
            vec![2.227563446120385],
            vec![2.6246054620212718],
            vec![2.447858975161667],
            vec![2.237357640550462],
            vec![2.0934313391366097],
            vec![3.333420423841205],
            vec![2.686746365701651],
            vec![2.6133102476651033],
            vec![2.479348669116296],
            vec![2.963689583658498],
            vec![2.4187689552578386],
            vec![1.7402971198064727],
            vec![2.1859493885988095],
            vec![2.710127360768819],
            vec![2.1599089641115228],
            vec![2.1339070837908545],
            vec![1.7414057861501107],
            vec![2.7053884800822203],
            vec![2.0500390858755906],
            vec![1.0452712048794233],
            vec![1.4549001430234874],
            vec![2.370954870925898],
            vec![1.6517644661658555],
            vec![1.2278959755109826],
            vec![2.0520445246086316],
            vec![2.3391891902465907],
            vec![2.1059071947015133],
            vec![2.0940952606686207],
            vec![1.4336383727971627],
            vec![2.4349485232892585],
            vec![2.070326684251644],
            vec![1.6822609547423861],
            vec![1.1248495603534834],
            vec![2.2241889358469344],
            vec![1.9546677186759482],
            vec![0.9155388907133534],
            vec![1.4939677743124844],
            vec![1.8635020162560347],
            vec![1.4922342853320678],
            vec![0.8218686846178962],
            vec![0.9829008219758708],
            vec![1.6113909939051807],
            vec![1.4219218997163696],
            vec![0.8466168975581255],
            vec![1.0958518277683575],
            vec![1.757656235365733],
            vec![0.9777544944392111],
            vec![0.8346588763591263],
            vec![0.7830090618919947],
            vec![1.1344487047852785],
            vec![0.3399650317783923],
            vec![-0.007620667611160592],
            vec![0.4005001159190091],
            vec![1.0852734560400088],
            vec![0.25194663124156025],
            vec![0.12013558004965719],
            vec![0.038796103126600714],
            vec![1.1295339852896424],
            vec![0.1555759903412322],
            vec![-0.2680103623120031],
            vec![0.0346930761249463],
            vec![0.7892883054682143],
            vec![0.4518868126500978],
            vec![0.174781420604139],
            vec![0.06530027535552485],
            vec![1.1732151637178023],
            vec![0.7726258068196201],
            vec![0.5143644649072533],
            vec![0.5971340921023655],
            vec![1.68133498584108],
            vec![1.2185144011770643],
            vec![1.4470371137399578],
            vec![1.4496130040713815],
            vec![1.6969619334871209],
            vec![1.4056482299096158],
            vec![1.1640066895062084],
            vec![0.7649655632989143],
            vec![1.5267045476010226],
        ];
        let j = run(
            "state_space_smoothing",
            &rows,
            &[5., 4., 5., 0.95, 2500., 1e-8],
        )
        .unwrap();
        assert!(
            (number(&j, "sigma_squared") - 0.08108991773316335).abs() < 2e-6,
            "{j}"
        );
        let par = values(&j, "smoothing_parameters");
        let want = [
            0.4361894792014159,
            0.10082997931708015,
            3.0294566311832646e-11,
            0.8000000010321984,
        ];
        for i in 0..4 {
            assert!((par[i] - want[i]).abs() < 0.0004, "{par:?} != {want:?}");
        }
        assert!(j.contains("\"converged\":true"), "{j}");
        let mut args = vec![5., 0.95];
        args.extend(values(&j, "prediction_model"));
        let k = run("state_space_predict", &[], &args).unwrap();
        check(&values(&j, "forecast_means"), &values(&k, "forecast_means"));
        check(
            &values(&j, "forecast_observation_variances"),
            &values(&k, "forecast_observation_variances"),
        );
    }
}
#[test]
fn ets_ann_forecast_variance_uses_shared_innovations_and_failed_optimizer_has_no_ic() {
    let rows = vec![
        vec![2.0000342874782],
        vec![1.8167048890869033],
        vec![2.166233669768416],
        vec![2.4319014093536797],
        vec![1.9145671885706752],
        vec![2.2319058081630043],
        vec![1.8208108446820492],
        vec![1.9874718189737965],
        vec![2.2978627440436132],
        vec![2.171892731841144],
        vec![1.5485318747658159],
        vec![2.371961926059055],
        vec![1.5798888805582152],
        vec![1.5054699738486075],
        vec![1.8688839526119838],
        vec![1.8693906291991589],
        vec![1.6827784209059515],
        vec![1.6436353611943137],
        vec![1.4872856545940125],
        vec![1.7898123463522735],
        vec![1.7817886754290673],
        vec![1.7298807419968385],
        vec![1.3939250911178402],
        vec![1.9804510339191803],
        vec![1.718139940995836],
        vec![2.00374632849518],
        vec![1.9397600805686943],
        vec![1.9440245159758172],
        vec![1.8919822230680545],
        vec![1.983588592247076],
        vec![1.8695526314740125],
        vec![1.8452449541305507],
        vec![2.27824081564765],
        vec![1.4402425369556728],
        vec![2.089198068409359],
        vec![1.7859874498349018],
        vec![1.84274268915339],
        vec![1.8408527381652366],
        vec![2.4155260626863333],
        vec![1.9704541634733614],
        vec![2.0758589539052767],
        vec![1.9784078013597042],
        vec![1.7903550166297706],
        vec![2.4187116629822465],
        vec![2.5477971654466813],
        vec![2.5818108683295393],
        vec![2.2634909206544567],
        vec![2.864083210376393],
        vec![2.1230315818103582],
        vec![2.503111216236104],
        vec![2.707902295132613],
        vec![2.380642221081157],
        vec![2.6518033667795993],
        vec![2.3657675398169404],
        vec![2.2822055282847358],
        vec![2.6210995515592215],
        vec![2.8963978093107046],
        vec![2.442351394553133],
        vec![2.9629388974197512],
        vec![2.634222261380734],
        vec![2.944554092838275],
        vec![2.2778255027858134],
        vec![2.2662486809267115],
        vec![2.8731271073756504],
    ];
    let j = run(
        "state_space_smoothing",
        &rows,
        &[0., 0., 4., 0.95, 2000., 1e-8],
    )
    .unwrap();
    let a = values(&j, "smoothing_parameters")[0];
    let s = number(&j, "sigma_squared");
    let v = values(&j, "forecast_observation_variances");
    for (i, value) in v.iter().enumerate() {
        close(*value, s * (1. + i as f64 * a * a));
    }
    let k = run(
        "state_space_smoothing",
        &rows,
        &[0., 0., 4., 0.95, 1., 1e-8],
    )
    .unwrap();
    assert!(k.contains("\"converged\":false"), "{k}");
    assert!(k.contains("\"aic\":null"), "{k}");
    let moved = rows
        .iter()
        .map(|r| vec![1e9 + 2. * r[0]])
        .collect::<Vec<_>>();
    let k = run(
        "state_space_smoothing",
        &moved,
        &[0., 0., 4., 0.95, 2000., 1e-8],
    )
    .unwrap();
    for (x, y) in values(&j, "forecast_means")
        .iter()
        .zip(values(&k, "forecast_means"))
    {
        assert!((y - (1e9 + 2. * x)).abs() < 1e-5);
    }
    assert!((number(&k, "sigma_squared") - 4. * s).abs() < 1e-6);
}

#[test]
fn portable_ets_scales_variance_without_an_overflowing_intermediate_and_rejects_underflow() {
    // k=1: level unchanged without future innovations, innovations gain alpha=.5.
    let mut model = vec![1., 1., 1., 0., 1e200, 1e-200, 0., 1., 1., 0.5];
    let mut args = vec![2., 0.95];
    args.extend(&model);
    let j = run("state_space_predict", &[], &args).unwrap();
    let v = values(&j, "forecast_observation_variances");
    assert!((v[0] / 1e200 - 1.).abs() < 1e-12);
    assert!((v[1] / 1e200 - 1.25).abs() < 1e-12);
    model[4] = 1e-300;
    model[5] = 1e-300;
    let mut args = vec![2., 0.95];
    args.extend(&model);
    assert!(run("state_space_predict", &[], &args).is_err());
}

#[test]
fn portable_forecast_accepts_one_ignored_scalar_row_for_the_native_nonempty_matrix_contract() {
    let args = [2., 0.95, 1., 0., 1., 1., 0.25, 1., 0.2, 1., 0.5];
    let j = run("state_space_predict", &[vec![0.]], &args).unwrap();
    check(&values(&j, "forecast_means"), &[1., 1.]);
    check(&values(&j, "forecast_observation_variances"), &[0.95, 1.15]);
    let k = run("state_space_predict", &[vec![1e200]], &args).unwrap();
    check(&values(&j, "forecast_means"), &values(&k, "forecast_means"));
    assert!(run("state_space_predict", &[vec![0.], vec![0.]], &args).is_err());
    assert!(run("state_space_predict", &[vec![f64::NAN]], &args).is_err());
}

fn two_state_external_covariance_args() -> Vec<f64> {
    vec![
        2., 1., 0.95, 0., 0., 1., 0., 0., 1., 1., 0., 0., 1., 1., 0., 0., 1., 1., 0., 1.,
    ]
}
#[test]
fn external_covariance_rejects_roundoff_asymmetry_instead_of_using_one_triangle() {
    for at in [6, 14] {
        let mut a = two_state_external_covariance_args();
        a[at] = 1e-13;
        assert!(run("kalman_smoother", &[vec![0., 1.]], &a).is_err());
    }
}
#[test]
fn external_process_covariance_rejects_small_negative_pivot() {
    let mut a = two_state_external_covariance_args();
    a[13] = -1e-13;
    assert!(run("kalman_smoother", &[vec![0., 1.]], &a).is_err());
}
#[test]
fn saved_process_covariance_cannot_hide_negative_latent_variance_behind_large_measurement_noise() {
    let mut v = two_state_external_covariance_args();
    v[..3].copy_from_slice(&[1., 0., 2.]);
    v[5] = 0.;
    v[13] = -1e-13;
    v[19] = 1e20;
    let mut a = vec![1., 0.95];
    a.extend(v);
    assert!(run("state_space_predict", &[vec![0.]], &a).is_err());
}
