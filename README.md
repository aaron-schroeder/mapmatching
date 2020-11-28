# mapmatching

> Intelligently matches a series of recorded GPS coordinates to a known network of trails.

[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![License](http://img.shields.io/:license-mit-blue.svg)](http://badges.mit-license.org)

<!--
[![Build Status](http://img.shields.io/travis/badges/badgerbadgerbadger.svg?style=flat-square)](https://travis-ci.org/badges/badgerbadgerbadger) [![Dependency Status](http://img.shields.io/gemnasium/badges/badgerbadgerbadger.svg?style=flat-square)](https://gemnasium.com/badges/badgerbadgerbadger) [![Coverage Status](http://img.shields.io/coveralls/badges/badgerbadgerbadger.svg?style=flat-square)](https://coveralls.io/r/badges/badgerbadgerbadger) [![Badges](http://img.shields.io/:badges-9/9-ff6799.svg?style=flat-square)](https://github.com/badges/badgerbadgerbadger)
-->


---

## Table of Contents

- [Motivation](#motivation)
- [The Map Matching Algorithm](#the-map-matching-algorithm)
- [Example](#example)
- [Project Status](#project-status)
- [Inspiration](#inspiration)
- [Contact](#contact)
- [License](#license)

---

## Motivation

GPS measurements suffer from location uncertainty. [The US Government performance standard for GPS](https://www.gps.gov/technical/ps/2008-SPS-performance-standard.pdf) indicated an average horizontal accuracy of 4 meters as of 2008. It is an amazing technological feat that a person holding a GPS-enabled phone or watch can be located so precisely by a system of satellites orbiting 20 million meters overhead. 

Unfortunately, certain applications demand a higher level of accuracy. Trail runners are keenly interested in the slope of the trails where their workouts take place. Pace and slope are by far the two most important factors that determine workout intensity. A 10-minute mile may not sound particularly intense, until considering that the runner performed this feat on an average 20 percent grade. At each point along a trail runner's GPS trace, accurate elevation data is key to understanding the slope of the trail, and by proxy, the runner's intensity. 

GPS fitness trackers do not directly collect elevation data. Instead, GPS coordinates are typically collected every second and input into a digital elevation model (DEM), which represents a simplified view of the earth surface elevation above sea level. Even with a perfect DEM, inaccurate GPS locations can generate a false elevation profile for a workout. 

This is less of a problem where the Earth's topography varies slowly and predictably. Think of gently rolling bluffs in Kansas; a few meter's horizontal difference would barely affect the elevation coordinate. But now, imagine a trail runner traversing a steep hillside - running across it, not up or down. In reality, this runner is likely traveling on a trail with a relatively predictable and shallow slope. But the GPS device doesn't know about the trail, and it may indicate that the runner is rapidly traveling uphill or downhill as inaccurate location measurements land on either side of the actual trail. Any attempt to infer the intensity of the runner's workout must correct for this obvious misrepresentation of the runner's work. 

If we are only interested in long averages of intensity, on the order of minutes, the average slope would not suffer too much from GPS inaccuracy; the distance used to calculate average slope would be much larger than the error. But what about 5-second intensity? This is where high-quality location data is vital. Until exercise scientists develop the instrumentation to directly measure runners' instantaneous power output, we must use a proxy - speed and slope.

So how do we improve the elevation data we are using? By using known trail locations. Let's assume we have a perfect DEM, one which completely realistically represents the Earth's surface. And we have GPS data collected on a trail run at one-second intervals. This project implements an algorithm to constrain an inaccurate GPS trace to a known trail network. The  methodology is described below.
<!--
In Boulder, where I live, I can get ahold of high-quality trail location coordinates that are the result of actual surveys. Not everyone has access to such accurate and granular data, but resources like the [Trail Run Project](https://www.trailrunproject.com/) provide GPX files with coordinates for many popular trails. 
-->

---

## The Map Matching Algorithm

This project implements a map matching algorithm using a Hidden Markov Model (HMM) described in [a 2009 paper by Paul Newson and John Krumm](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.187.5145&rep=rep1&type=pdf). The linked paper describes the algorithm in detail. 

In a nutshell, the algorithm takes a GPS trace and finds the most likely corresponding path through a trail network. A path is considered "more likely" if each GPS point is close to the corresponding matched point, and if the point-to-point distance along the matched route is similar to the distance along the GPS route. In other words, the algorithm uses the GPS data as a starting point and attempts to produce a similar route on the trail grid. Since GPS errors are not too large, this often means that a given GPS point is matched to the closest point on the trail grid. But often, trails are winding and may even double back on themselves, meaning a simple closest-point map matching algorithm could jump forward or backward on the route by an unrealistic amount from one second to the next. The algorithm resolves this problem by getting a little smart and considering how far the user moved between GPS points.

---

## Example
```python
import networkx as nx
from geopy.distance import distance.great_circle as distance
from mapmatching import mapmatch

# Build the trail network graph using networkx
latlons_trailnet = [[-105.000, 40.000], [-105.000, 40.100], [-105.200, 40.100]]
graph_trailnet = nx.Graph()
for i in range(0, len(latlons_trailnet)):
    graph.add_node(
        i,
        pos=(latlons_trailnet[i][1], latlons_trailnet[i][0])
    )
for i in range(0, len(latlons_trailnet)-1):
    graph.add_edge(
        i,
        i+1,
        key=i,
        length=distance(
            latlons_trailnet[i][::-1],
            latlons_trailnet[i+1][::-1]).meters
    )

# Build gps lists/arrays (typically these would be loaded from a file)
latlons_gps = [[-105.000, 40.00000], [-105.000, 40.00003], 
               [-105.000, 40.00006], [-105.000, 40.00009],
               [-105.000, 40.00012], [-105.000, 40.00015],
               [-105.000, 40.00018], [-105.000, 40.00021],
               [-105.000, 40.00024], [-105.000, 40.00024]]
times_gps = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]               

# Feed the input into the map matching algorithm
points_matched, emission_points_matched = mapmatch(
    latlons_gps,
    times_gps,
    graph_trailnet)
```

---

## Project Status

### Complete

- Create Python package. Make pretty input-output later.

### Current Activities

- Implement a series of tests to ensure functionality as development progresses.

- Streamline input so user can be more hands-off.

#### Benchmarking and Optimization

- Benchmark algorithm performance (speed and accuracy). Generate fake user GPS coordinates by sampling along a route and simulating GPS drift. Document algorithm behavior under various conditions by varying: 
   - simulated GPS sampling rate
   - magnitude of simulated GPS error (standard deviation)
   - trail network segment lengths relative to GPS segment lengths
   - algorithm parameters: SIGMA_Z, BETA, DIST_ERR_TOL, etc

### Future Work

- Enhance the speed of the Viterbi algorithm by minimizing the computational load as much as possible without degrading performance. This requires the results of the performance tests described above.

- Compare algorithm performance to publicly available map matching services.
   - [Mapbox](https://docs.mapbox.com/api/navigation/#map-matching)
   - Others?

- Implement a version of [Dijkstra's algorithm](https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm), described in the [documentation for another map matching project](https://github.com/valhalla/valhalla/blob/master/docs/meili/algorithms.md) by Valhalla, as an alternate to the Viterbi algorithm which is currently in use.

- Compare the speed and accuracy of the two algorithm options.

- Integrate OpenStreetMap to provide trail network data.

- Integrate with my website's route database.

- Create a page on my website that allows users to match their own gps coordinates.

---

<!--
## FAQ
- **How do I do *specifically* so and so?**
    - No problem! Just do this.
-->
<!--
## People

| <a href="https://github.com/aaron-schroeder" target="_blank">**Aaron Schroeder**</a> |
| :---: |
| [![Aaron Schroeder](https://avatars0.githubusercontent.com/u/39806580?v=4&s=200)](https://github.com/aaron-schroeder) |
| <a href="https://github.com/aaron-schroeder" target="_blank">`github.com/aaron-schroeder`</a> |

---
-->

## Inspiration

- [A Developer Diary](http://www.adeveloperdiary.com/data-science/machine-learning/implement-viterbi-algorithm-in-hidden-markov-model-using-python-and-r/) for helping understand the nuts and bolts of the Viterbi algorithm in Python.

- [Mapzen](https://www.mapzen.com/blog/map-matching-validation/), original creators of [Valhalla](https://www.interline.io/valhalla/), for ideas for benchmarking and future work.

- [Graphhopper](https://github.com/graphhopper/map-matching) for description of algorithm implementation.

---

## Contact

Reach out to me at one of the following places!

- Website: <a href="https://trailzealot.com" target="_blank">trailzealot.com</a>
- LinkedIn: <a href="https://www.linkedin.com/in/aarondschroeder/" target="_blank">linkedin.com/in/aarondschroeder</a>
- Twitter: <a href="https://twitter.com/trailzealot" target="_blank">@trailzealot</a>
- Instagram: <a href="https://instagram.com/trailzealot" target="_blank">@trailzealot</a>
- GitHub: <a href="https://github.com/aaron-schroeder" target="_blank">github.com/aaron-schroeder</a>

---

<!--
## Support

[![Buy me a coffee](https://bmc-cdn.nyc3.digitaloceanspaces.com/BMC-button-images/custom_images/orange_img.png)](https://www.buymeacoffee.com/OpPEjlGFC)

[![Support via Patreon](https://c5.patreon.com/external/logo/become_a_patron_button@2x.png)](https://www.patreon.com/amitmerchant)

[![Support via Gratipay](https://cdn.rawgit.com/gratipay/gratipay-badge/2.3.0/dist/gratipay.png)](https://gratipay.com/fvcproductions/)

---
-->

## License

[![License](http://img.shields.io/:license-mit-blue.svg)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**
- Copyright 2019 © <a href="https://trailzealot.com/about" target="_blank">Aaron Schroeder</a>.

<!--
## Installation

- All the `code` required to get started
- Images of what it should look like

### Clone

- Clone this repo to your local machine using `https://github.com/EricSchraider/mapmatching`

### Setup

- If you want more syntax highlighting, format your code like this:

> update and install this package first

```shell
$ brew update
$ brew install fvcproductions
```

> now install npm and bower packages

```shell
$ npm install
$ bower install
```

---
-->

<!--
## Contributing

> To get started...

### Step 1

- **Option 1**
    - 🍴 Fork this repo!

- **Option 2**
    - 👯 Clone this repo to your local machine using `https://github.com/joanaz/HireDot2.git`

### Step 2

- **HACK AWAY!** 🔨🔨🔨

### Step 3

- 🔃 Create a new pull request using <a href="https://github.com/joanaz/HireDot2/compare/" target="_blank">`https://github.com/joanaz/HireDot2/compare/`</a>.

---
-->


<!--
## Features
## Documentation (Optional)
## Tests (Optional)
-->
