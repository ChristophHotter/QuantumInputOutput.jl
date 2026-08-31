window.BENCHMARK_DATA = {
  "lastUpdate": 1787856894764,
  "repoUrl": "https://github.com/qojulia/QuantumInputOutput.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b0a8e8b217daa7690878219854dfc71a81548354",
          "message": "Merge pull request #5 from ChristophHotter/static\n\nrefactor: More type-stable and concrete implementation",
          "timestamp": "2026-04-13T16:30:58-04:00",
          "tree_id": "ae0b973479310ec1f05e246528ce8fa6d0028e04",
          "url": "https://github.com/ChristophHotter/QuantumInputOutput.jl/commit/b0a8e8b217daa7690878219854dfc71a81548354"
        },
        "date": 1776113014126,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 288322940,
            "unit": "ns",
            "extra": "gctime=44054876.5\nmemory=373256856\nallocs=2779246\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 20791,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 143355,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47880\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 460.99,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":200,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 2893.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 567725595,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53142384\nallocs=1452354\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 449433,
            "unit": "ns",
            "extra": "gctime=0\nmemory=986320\nallocs=554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 505557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141568\nallocs=451\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 363835,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 144.1439075630252,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":952,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 405367,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 138.12225705329155,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":957,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 87811.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=129232\nallocs=291\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 10378.0625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 67601,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84896\nallocs=377\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 100471,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26288\nallocs=724\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 32929,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9280\nallocs=257\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 33711,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9312\nallocs=297\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 20250,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106680\nallocs=30\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 15423,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85648\nallocs=20\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 444615.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45648\nallocs=1049\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 3275953,
            "unit": "ns",
            "extra": "gctime=0\nmemory=648920\nallocs=8471\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 799647,
            "unit": "ns",
            "extra": "gctime=0\nmemory=76624\nallocs=1856\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e8e1d052b100618271bf2f24916b7ddeb7ca2aa5",
          "message": "Merge pull request #6 from ChristophHotter/code_quality\n\nCode quality",
          "timestamp": "2026-04-16T11:09:11-04:00",
          "tree_id": "80c04f595dba2c605dc7ee58507b7a77c17ad308",
          "url": "https://github.com/ChristophHotter/QuantumInputOutput.jl/commit/e8e1d052b100618271bf2f24916b7ddeb7ca2aa5"
        },
        "date": 1776352779330,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 206796793,
            "unit": "ns",
            "extra": "gctime=31767955\nmemory=373256856\nallocs=2779246\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 18975,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 132827,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47880\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 380.2682926829268,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":205,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 2689.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 588152097,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53129856\nallocs=1452109\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 440810.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=986320\nallocs=554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 490457.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141568\nallocs=451\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 380247,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 126.79769392033543,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":954,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 416925.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 133.66041666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":960,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 77635,
            "unit": "ns",
            "extra": "gctime=0\nmemory=129232\nallocs=291\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 8684.3125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 62817,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84896\nallocs=377\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 126656,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26288\nallocs=724\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 39593,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9280\nallocs=257\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 35916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9312\nallocs=297\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 21384.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106680\nallocs=30\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 15960,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85648\nallocs=20\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 501523.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45648\nallocs=1049\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 3560007.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=642040\nallocs=8403\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 892202,
            "unit": "ns",
            "extra": "gctime=0\nmemory=77248\nallocs=1863\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e8824b7fb6c61373365bd1a8add7567a8cca722f",
          "message": "Merge pull request #9 from qojulia/new_tests\n\nadded more tests",
          "timestamp": "2026-04-17T14:22:48-04:00",
          "tree_id": "df9ef2efa2af16812f72849099cbe41f0d6e4641",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/e8824b7fb6c61373365bd1a8add7567a8cca722f"
        },
        "date": 1776450504249,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 256011981,
            "unit": "ns",
            "extra": "gctime=39761984.5\nmemory=373256856\nallocs=2779246\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 19016,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 140923,
            "unit": "ns",
            "extra": "gctime=0\nmemory=47880\nallocs=547\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 370.15533980582524,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":206,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 2715,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 605594816,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53129856\nallocs=1452109\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 456150,
            "unit": "ns",
            "extra": "gctime=0\nmemory=986320\nallocs=554\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 510953,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141568\nallocs=451\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 387638.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 137.75729166666667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":960,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 419272.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 154.90688912809472,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":929,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 78446,
            "unit": "ns",
            "extra": "gctime=0\nmemory=129232\nallocs=291\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 9179.625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 62707,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84896\nallocs=377\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 130644,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26288\nallocs=724\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 41698,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9280\nallocs=257\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 36017,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9312\nallocs=297\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 20779,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106680\nallocs=30\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 15950,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85648\nallocs=20\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 518337,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45648\nallocs=1049\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 3930011.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=642040\nallocs=8403\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 904327,
            "unit": "ns",
            "extra": "gctime=0\nmemory=77248\nallocs=1863\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2a155abd41e1239454750ab1e2796b75643caf1b",
          "message": "Merge pull request #10 from qojulia/adjoint_ops_translate\n\nAdjoint ops translate",
          "timestamp": "2026-05-31T18:01:33-04:00",
          "tree_id": "9460da8d2f545440df716749bbb816118cdd52d8",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/2a155abd41e1239454750ab1e2796b75643caf1b"
        },
        "date": 1780265649484,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 245918285,
            "unit": "ns",
            "extra": "gctime=34038238\nmemory=373176792\nallocs=2777826\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 17328,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 121594.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48024\nallocs=548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 311.7795275590551,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":254,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 2145.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 531590721,
            "unit": "ns",
            "extra": "gctime=0\nmemory=53138208\nallocs=1452247\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 437219,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984640\nallocs=505\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 510509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141296\nallocs=438\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 380661,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 127.10699373695198,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":958,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 404046.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 130.90388655462186,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":952,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 65417.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=129232\nallocs=291\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 7646.375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 52371.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=84896\nallocs=377\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 98119.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=26288\nallocs=724\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 31346,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9280\nallocs=257\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 32864,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9312\nallocs=297\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 23670,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106680\nallocs=30\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 18372,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85648\nallocs=20\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 498433.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=45648\nallocs=1049\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 3576277,
            "unit": "ns",
            "extra": "gctime=0\nmemory=664256\nallocs=8329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 893640,
            "unit": "ns",
            "extra": "gctime=0\nmemory=76096\nallocs=1866\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "009b4390f1c525338cc317c4aa00bf2b94d759de",
          "message": "Merge pull request #15 from qojulia/sqav0.6\n\nmove to SQA v0.8",
          "timestamp": "2026-06-25T15:07:13+02:00",
          "tree_id": "b0831363b4bddb095636449ccd4534840f898539",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/009b4390f1c525338cc317c4aa00bf2b94d759de"
        },
        "date": 1782393779581,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 282088691,
            "unit": "ns",
            "extra": "gctime=34289421\nmemory=373176344\nallocs=2777821\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13871,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 73451,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48024\nallocs=548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 109.43257820927724,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":927,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 220.56716417910448,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":469,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 78759,
            "unit": "ns",
            "extra": "gctime=0\nmemory=74016\nallocs=1243\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 425233,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984640\nallocs=505\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 490913,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141296\nallocs=438\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 353034,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 91.06792058516197,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":957,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 389409,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 90.9633891213389,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":956,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 37086,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3542.875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 21893,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54064\nallocs=287\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 8153,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23152\nallocs=241\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 3043.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10480\nallocs=100\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 7023.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20704\nallocs=182\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 19449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106680\nallocs=30\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 14712,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85648\nallocs=20\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 32920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20928\nallocs=396\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 502680,
            "unit": "ns",
            "extra": "gctime=0\nmemory=471600\nallocs=3647\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 86711,
            "unit": "ns",
            "extra": "gctime=0\nmemory=40336\nallocs=914\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ac622d7faefae370aa0e53876a202e66702a336e",
          "message": "Bump version from 0.1.0 to 0.2.0 (#26)",
          "timestamp": "2026-06-30T15:45:44-04:00",
          "tree_id": "7aebf46a253f1bd4794ad9f452d241d1267b4e6c",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/ac622d7faefae370aa0e53876a202e66702a336e"
        },
        "date": 1782849583216,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 246217008,
            "unit": "ns",
            "extra": "gctime=26568491\nmemory=373176344\nallocs=2777821\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13305,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 72195,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48024\nallocs=548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 96.31296101159114,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":949,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 211.65142857142857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":525,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 81924,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73904\nallocs=1241\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 416737,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984640\nallocs=505\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 460660,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141296\nallocs=438\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 366424,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 88.25833333333334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":960,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 396551,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 92.27882599580713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":954,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 33602,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3311.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 21169,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54064\nallocs=287\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7818,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23152\nallocs=241\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2690.6666666666665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10480\nallocs=100\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6410,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20704\nallocs=182\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 19055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106680\nallocs=30\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 14527,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85648\nallocs=20\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 38251,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20928\nallocs=396\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 569613,
            "unit": "ns",
            "extra": "gctime=0\nmemory=504152\nallocs=3657\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 114223,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39696\nallocs=907\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "39d5b8b12dd337cfbf45d48d25ac5c6400ca63ec",
          "message": "Super example (#28)",
          "timestamp": "2026-07-09T12:09:14-04:00",
          "tree_id": "ae75e3f5ce12c7b1653bd0b2376913df5e29c8c7",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/39d5b8b12dd337cfbf45d48d25ac5c6400ca63ec"
        },
        "date": 1783614133093,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 246696455,
            "unit": "ns",
            "extra": "gctime=26234324\nmemory=373176344\nallocs=2777821\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13045,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 72716,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48024\nallocs=548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 99.35911016949153,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":944,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 201.46101694915254,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":590,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 83545,
            "unit": "ns",
            "extra": "gctime=0\nmemory=73904\nallocs=1241\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 407841,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984640\nallocs=505\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 456171,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141296\nallocs=438\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 364751,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 103.44893617021276,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":940,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 397302,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 92.09789473684211,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":950,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 33863,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3214.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 21911,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54064\nallocs=287\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7968,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23152\nallocs=241\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2811.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10480\nallocs=100\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6492.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20704\nallocs=182\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 19115,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106680\nallocs=30\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 14637,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85648\nallocs=20\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 37991,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20928\nallocs=396\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 573420,
            "unit": "ns",
            "extra": "gctime=0\nmemory=471600\nallocs=3647\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 115035,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39696\nallocs=907\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ccbe25777c92d0bfe5ad4fe2b214a0442f34b73d",
          "message": "feat: update to SQA v0.9 (#30)",
          "timestamp": "2026-07-12T18:12:28+02:00",
          "tree_id": "d3b7f38ec6f53f65f8a9e9745244c9dc31127ddd",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/ccbe25777c92d0bfe5ad4fe2b214a0442f34b73d"
        },
        "date": 1783873551604,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 246747052,
            "unit": "ns",
            "extra": "gctime=26984577\nmemory=373176344\nallocs=2777821\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13044,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 72365,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48024\nallocs=548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 96.00429184549357,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":932,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 212.87857142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":560,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 18214,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41872\nallocs=562\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 411919,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984640\nallocs=505\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 459147,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141296\nallocs=438\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 364340,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 95.79284963196635,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":951,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 397021,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 91.01567398119123,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":957,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 34224,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3402.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 21590,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54064\nallocs=287\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7955,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23152\nallocs=241\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2824.222222222222,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10480\nallocs=100\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6395.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20704\nallocs=182\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 19326,
            "unit": "ns",
            "extra": "gctime=0\nmemory=106680\nallocs=30\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 14637,
            "unit": "ns",
            "extra": "gctime=0\nmemory=85648\nallocs=20\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 41988,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22800\nallocs=437\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 580312,
            "unit": "ns",
            "extra": "gctime=0\nmemory=475776\nallocs=3737\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 119603,
            "unit": "ns",
            "extra": "gctime=0\nmemory=42464\nallocs=958\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "8e7b370fb77f5b7073a7d709a9ae5c2610044059",
          "message": "Bump version from 0.3.0 to 0.4.0 (#33)",
          "timestamp": "2026-07-23T18:50:54+02:00",
          "tree_id": "ded347cab05debc33d3d4fd555ba5db7adec866a",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/8e7b370fb77f5b7073a7d709a9ae5c2610044059"
        },
        "date": 1784826351299,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 162497494,
            "unit": "ns",
            "extra": "gctime=16810722\nmemory=197626952\nallocs=576686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 10527,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 60411,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48024\nallocs=548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 86.58896982310094,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":961,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 169.01066666666668,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":750,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 14212,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 292091,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984560\nallocs=504\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 353395,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141216\nallocs=437\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 231499,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 62.66836734693877,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":980,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 262018,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 62.843877551020405,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":980,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 28790,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 2866.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 16303,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54064\nallocs=287\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 6131.833333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22416\nallocs=221\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2259.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10192\nallocs=92\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20704\nallocs=182\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 5744.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 5709.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 25245,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21728\nallocs=410\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 425183,
            "unit": "ns",
            "extra": "gctime=0\nmemory=418944\nallocs=3729\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 71412,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41456\nallocs=949\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6a468fd227b4ea6acfb105d5aa5460704eb0f5e6",
          "message": "build: move to SQA v0.10 (#32)",
          "timestamp": "2026-07-23T18:50:21+02:00",
          "tree_id": "5d80040d31aeffe48b8a3e7dd10833a7810b8ba2",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/6a468fd227b4ea6acfb105d5aa5460704eb0f5e6"
        },
        "date": 1784826407138,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 183117800,
            "unit": "ns",
            "extra": "gctime=13662582\nmemory=197626952\nallocs=576686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13104,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 73197,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48024\nallocs=548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 96.30126849894292,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":946,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 204.60943396226415,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":530,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 17773,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 414315,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984560\nallocs=504\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 464319,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141216\nallocs=437\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 369711,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 89.57306889352819,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":958,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 403004,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 90.81799163179916,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":956,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 34375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3233.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 21570,
            "unit": "ns",
            "extra": "gctime=0\nmemory=54064\nallocs=287\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7576.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22416\nallocs=221\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2702.777777777778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10192\nallocs=92\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6590.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20704\nallocs=182\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 8616,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 9077,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 39484,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21728\nallocs=410\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 596616,
            "unit": "ns",
            "extra": "gctime=0\nmemory=418944\nallocs=3729\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 122599,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41456\nallocs=949\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6656d37151f8a4113af549324c2145c72931b19c",
          "message": "fix: make SLH scattering field concretely typed (#34)",
          "timestamp": "2026-08-06T10:01:04+02:00",
          "tree_id": "8b427c72a02491dc0b43a173c2f92d00696c7ab9",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/6656d37151f8a4113af549324c2145c72931b19c"
        },
        "date": 1786004226933,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 180346610,
            "unit": "ns",
            "extra": "gctime=11829440\nmemory=197626952\nallocs=576686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13085,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16560\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 72897,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48024\nallocs=548\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 98.28224101479915,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":946,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 203.49734513274336,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":565,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 17252,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 407135,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984560\nallocs=504\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 455866,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141216\nallocs=437\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 364786,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561320\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 88.73695198329854,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":958,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 396505,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721456\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 95.5099894847529,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":951,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 34244,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3226,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 30828,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 6925,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2157.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 5953.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 9067,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 9157,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 39404,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21728\nallocs=410\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 589577,
            "unit": "ns",
            "extra": "gctime=0\nmemory=418848\nallocs=3726\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 125606,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41424\nallocs=948\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a1d1e7e5253dd99158d332f845d223116fcbb5a1",
          "message": "build(deps): bump ODE (#35)",
          "timestamp": "2026-08-06T11:40:16+02:00",
          "tree_id": "95f6cdf2bdc445ffb547d28b8d600cc833fd7041",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/a1d1e7e5253dd99158d332f845d223116fcbb5a1"
        },
        "date": 1786009909952,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 151584403,
            "unit": "ns",
            "extra": "gctime=14776593\nmemory=198187032\nallocs=580186\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 10737,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 55674,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49032\nallocs=556\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 92.19392033542977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":954,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 186.8048780487805,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":656,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 15153,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 282795,
            "unit": "ns",
            "extra": "gctime=0\nmemory=985664\nallocs=512\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 378339,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1142144\nallocs=445\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 262444,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 85.847022587269,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":974,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 303837,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 91.02055498458377,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":973,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 30546,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 2977,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 25047,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 6457.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 1949.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 5595,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 8095.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 8215.666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 26760,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21728\nallocs=410\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 426942,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419968\nallocs=3726\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 81572,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41424\nallocs=948\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "dab446f829193b8c444767a7f57beab9bb1a4dca",
          "message": "fix: coupling_matrix to support complex couplings (#40)",
          "timestamp": "2026-08-06T13:17:14+02:00",
          "tree_id": "071999f9d7c8fd8477d20ddefd4abb256ceab070",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/dab446f829193b8c444767a7f57beab9bb1a4dca"
        },
        "date": 1786015856051,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 185359501,
            "unit": "ns",
            "extra": "gctime=13449344\nmemory=198187032\nallocs=580186\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13324,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 68980,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49032\nallocs=556\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 96.26027397260275,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":949,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 194.76426426426426,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":666,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 17653,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 415828,
            "unit": "ns",
            "extra": "gctime=0\nmemory=985664\nallocs=512\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 475458,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1142144\nallocs=445\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 370696,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 89.67640918580376,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":958,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 404718,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 102.32161874334399,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":939,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 34905,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3286.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 30958,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7018,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2194.1111111111113,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6193.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 8436,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 8616,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 45235,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21728\nallocs=410\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 604071,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419968\nallocs=3726\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 131576,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41424\nallocs=948\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "76204328772c9699b9ad3bb8390731ce44efb36a",
          "message": "Pulse coupling (#43)",
          "timestamp": "2026-08-09T11:16:28-04:00",
          "tree_id": "f94bf16fdd806b23254e3ea21a5063ef18baa2d7",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/76204328772c9699b9ad3bb8390731ce44efb36a"
        },
        "date": 1786289467494,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 182943988,
            "unit": "ns",
            "extra": "gctime=13812646\nmemory=198187128\nallocs=580186\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 14176,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 66143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=49032\nallocs=556\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 96.2601880877743,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":957,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 191.75118858954042,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":631,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 17964,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 407237,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984912\nallocs=494\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 457893,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141392\nallocs=427\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 368086,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 91.47803347280335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":956,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 402400,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 90.99686192468619,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":956,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 34223,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3236.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 30337,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 6953,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2236.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 5947,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 9388,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 9387,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 40095,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21728\nallocs=410\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 587876,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419968\nallocs=3726\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 123259,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41424\nallocs=948\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4220857888bfd205ecbc1d56532697ccf55274d9",
          "message": "rename lindblad to jump_operator (#44)",
          "timestamp": "2026-08-20T19:03:33-04:00",
          "tree_id": "5b4e6cac08f3258e743951258c7df0f47f3393d1",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/4220857888bfd205ecbc1d56532697ccf55274d9"
        },
        "date": 1787268272346,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 185881229,
            "unit": "ns",
            "extra": "gctime=12007264\nmemory=198107128\nallocs=578686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13545,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 69240,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48872\nallocs=553\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 96.90736842105264,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":950,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 190.49167927382754,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":661,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 18465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 415657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984768\nallocs=491\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 469136,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141248\nallocs=424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 370152,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 98.80147835269271,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":947,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 401781,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 96.37327677624603,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":943,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 34825,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3432.625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 30136,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7161,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2187.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6887.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 9318,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 9618,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 55344,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21632\nallocs=404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 713092,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419488\nallocs=3664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 162073,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41488\nallocs=940\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "eda9bd34cdbd6a618385a5efc8e7dbd8f802fb8c",
          "message": "Docu issue (#45)",
          "timestamp": "2026-08-21T09:23:51-04:00",
          "tree_id": "63b7ac6d0283f08b74d12deb07a6144fd1be2433",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/eda9bd34cdbd6a618385a5efc8e7dbd8f802fb8c"
        },
        "date": 1787319221272,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 161952373,
            "unit": "ns",
            "extra": "gctime=13908010\nmemory=198107128\nallocs=578686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 10796,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 52509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48872\nallocs=553\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 83.37927461139897,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":965,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 166.76923076923077,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":767,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 15083,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 290108,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984768\nallocs=491\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 391521,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141248\nallocs=424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 272761,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 73.1726618705036,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":973,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 312130,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 73.17368961973278,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":973,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 29935,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 2958.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 25108,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 6179.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 1925.111111111111,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 5720.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":6,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 7321,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 7110.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 40941,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21632\nallocs=404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 503800,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419488\nallocs=3664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 117497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41488\nallocs=940\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "752e72b50bacff335202c08e3e8f926dbda23a34",
          "message": "Bump version from 0.4.1 to 0.5.0 (#46)",
          "timestamp": "2026-08-21T10:14:12-04:00",
          "tree_id": "f1b9df3f5a7ce33b8182264d0090b79542d70274",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/752e72b50bacff335202c08e3e8f926dbda23a34"
        },
        "date": 1787322339408,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 199888911,
            "unit": "ns",
            "extra": "gctime=11350574\nmemory=198107128\nallocs=578686\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13941,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 69555,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48872\nallocs=553\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 103.27234042553191,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":940,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 208.02522522522523,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":555,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 19118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 434726,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984768\nallocs=491\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 491151,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141248\nallocs=424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 353343,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 94.42752100840336,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":952,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 389989,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 98.36947368421053,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":950,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 37607,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3655.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 30506,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7501.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2383.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6745,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 9555,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 9724,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 56315,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21632\nallocs=404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 604792,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419488\nallocs=3664\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 153922,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41488\nallocs=940\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "0680161bea37096e4261872d9740dabd4299ccb2",
          "message": "change L to J (#48)",
          "timestamp": "2026-08-24T11:16:14-04:00",
          "tree_id": "d23772899f78212809e5460cb1d78ab0cf719d9f",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/0680161bea37096e4261872d9740dabd4299ccb2"
        },
        "date": 1787585438816,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 146163489,
            "unit": "ns",
            "extra": "gctime=12789225\nmemory=198243128\nallocs=582186\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 10766,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 51308,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48872\nallocs=553\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 83.38278008298755,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":964,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 160.26658163265307,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":784,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 15474,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 279350,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984768\nallocs=491\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 380213,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141248\nallocs=424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 264088,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 74.21707818930041,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":972,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 303497,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 73.67728674203494,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":973,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 31137,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide J(t)",
            "value": 2875.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 24337,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 6173.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 1904.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":10,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 5215.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 7454.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity J(t)",
            "value": 7188.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 25939,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21632\nallocs=404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+J",
            "value": 416186,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419296\nallocs=3658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 74902,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41424\nallocs=938\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "918993186964d9df895c4f518c2c49ff0f518645",
          "message": "G1 hermitian (#49)",
          "timestamp": "2026-08-24T12:16:59-04:00",
          "tree_id": "4ee75ed8950f6e11a392b6ba1bc9bb8a2d258fd4",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/918993186964d9df895c4f518c2c49ff0f518645"
        },
        "date": 1787589211027,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 182102247,
            "unit": "ns",
            "extra": "gctime=13640309\nmemory=198244440\nallocs=582189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 14081,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 65578,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48872\nallocs=553\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 107.31969860064585,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":929,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 211.15,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":540,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 18758,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 426151,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984768\nallocs=491\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 484689,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141248\nallocs=424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 353351,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 93.94432773109244,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":952,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 389676,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 93.81970649895179,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":954,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 36836,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide J(t)",
            "value": 3647.875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 30165,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7343.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2343.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6157.2,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 8122,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity J(t)",
            "value": 8255.666666666666,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 32108,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21632\nallocs=404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+J",
            "value": 515045,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419296\nallocs=3658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 92309,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41424\nallocs=938\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "04404cde16a2d66afbb3c1e19124473aa252c33d",
          "message": "Bump version from 0.5.0 to 0.5.1 (#50)",
          "timestamp": "2026-08-24T12:19:04-04:00",
          "tree_id": "5eb0a182a38e9f7049044e514285afa51947229a",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/04404cde16a2d66afbb3c1e19124473aa252c33d"
        },
        "date": 1787589361564,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 187581210,
            "unit": "ns",
            "extra": "gctime=14782681\nmemory=198244440\nallocs=582189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13920,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 65136,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48872\nallocs=553\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 107.31808510638298,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":940,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 251.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":550,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 18888,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 430796,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984768\nallocs=491\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 490004,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141248\nallocs=424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 355995,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 93.82563025210084,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":952,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 392359,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 94.4359243697479,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":952,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 37806,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide J(t)",
            "value": 3679.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 31026,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7781.333333333333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":3,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2456.8888888888887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6351.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 9494,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity J(t)",
            "value": 9574,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 32608,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21632\nallocs=404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+J",
            "value": 527168,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419296\nallocs=3658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 94510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41424\nallocs=938\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "235551c19b135065707567314884e26e5c5d731d",
          "message": "Bump version from 0.5.1 to 0.5.2 (#52)",
          "timestamp": "2026-08-24T16:09:51-04:00",
          "tree_id": "daf28f3bd07769ad3411c69fb9e0744b8d0dadf3",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/235551c19b135065707567314884e26e5c5d731d"
        },
        "date": 1787602568070,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 182372099,
            "unit": "ns",
            "extra": "gctime=13225580\nmemory=198244440\nallocs=582189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13435,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 66393,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48872\nallocs=553\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 96.13277133825079,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":949,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 187.40712074303406,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":646,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 18154,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 422305,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984768\nallocs=491\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 451238,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141248\nallocs=424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 369928,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 97.4968287526427,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":946,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 390054,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 96.97993664202745,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":947,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 36338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide J(t)",
            "value": 3516.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 31018,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 6713,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2315.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 5961.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 9096,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity J(t)",
            "value": 9137,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 39453,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21632\nallocs=404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+J",
            "value": 610655,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419296\nallocs=3658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 123821,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41424\nallocs=938\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "chr.hotter@gmail.com",
            "name": "Christoph Hotter",
            "username": "ChristophHotter"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1e3ffbe6187a4af40db2d8f67e4786a2f4719470",
          "message": "Revert jump operator rename (#54)",
          "timestamp": "2026-08-27T14:35:36-04:00",
          "tree_id": "f46fc69486d3c524f3a7e554024d53c904fa13bf",
          "url": "https://github.com/qojulia/QuantumInputOutput.jl/commit/1e3ffbe6187a4af40db2d8f67e4786a2f4719470"
        },
        "date": 1787856891042,
        "tool": "julia",
        "benches": [
          {
            "name": "Correlations/two-time/single photon cavity",
            "value": 186364343,
            "unit": "ns",
            "extra": "gctime=13992185\nmemory=198244440\nallocs=582189\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/analytical (2 equal modes)",
            "value": 13530,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16784\nallocs=58\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coefficient matrix M/numerical (ODE)",
            "value": 65238,
            "unit": "ns",
            "extra": "gctime=0\nmemory=48872\nallocs=553\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/2 modes",
            "value": 107.32021276595745,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":940,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/coupling matrix evaluation/4 modes",
            "value": 211.1375,
            "unit": "ns",
            "extra": "gctime=0\nmemory=0\nallocs=0\nparams={\"evals\":560,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Interaction Picture/operator substitution/TLS cascade",
            "value": 18708,
            "unit": "ns",
            "extra": "gctime=0\nmemory=39824\nallocs=498\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/input 2 modes",
            "value": 420742,
            "unit": "ns",
            "extra": "gctime=0\nmemory=984768\nallocs=491\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/multi-pulse/output 2 modes",
            "value": 478448,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1141248\nallocs=424\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input",
            "value": 350436,
            "unit": "ns",
            "extra": "gctime=0\nmemory=561352\nallocs=32\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/input Gaussian",
            "value": 94.01570680628272,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":955,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output",
            "value": 388142,
            "unit": "ns",
            "extra": "gctime=0\nmemory=721488\nallocs=35\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Pulse Couplings/single-pulse/output Gaussian",
            "value": 94.2565445026178,
            "unit": "ns",
            "extra": "gctime=0\nmemory=112\nallocs=3\nparams={\"evals\":955,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide H(t)",
            "value": 36655,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130320\nallocs=321\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/closure evaluation/2-QD waveguide L(t)",
            "value": 3585.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16168\nallocs=29\nparams={\"evals\":8,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/numeric/2-QD waveguide composition",
            "value": 30085,
            "unit": "ns",
            "extra": "gctime=0\nmemory=56704\nallocs=329\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/3-cavity cascade",
            "value": 7413.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22016\nallocs=205\nparams={\"evals\":4,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/concatenate + cascade",
            "value": 2381.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9952\nallocs=82\nparams={\"evals\":9,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "SLH Algebra/symbolic/feedback OPO loop",
            "value": 6315.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20368\nallocs=168\nparams={\"evals\":5,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity H(t)",
            "value": 9103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=59344\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/closure evaluation/3-cavity L(t)",
            "value": 9204,
            "unit": "ns",
            "extra": "gctime=0\nmemory=63824\nallocs=55\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/static/atom-cavity",
            "value": 32189,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21632\nallocs=404\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/3-cavity H+L",
            "value": 518008,
            "unit": "ns",
            "extra": "gctime=0\nmemory=419296\nallocs=3658\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          },
          {
            "name": "Translation/time-dependent/atom-cavity",
            "value": 91868,
            "unit": "ns",
            "extra": "gctime=0\nmemory=41424\nallocs=938\nparams={\"evals\":1,\"evals_set\":false,\"gcsample\":false,\"gctrial\":true,\"memory_tolerance\":0.01,\"overhead\":0,\"samples\":10000,\"seconds\":5,\"time_tolerance\":0.05}"
          }
        ]
      }
    ]
  }
}