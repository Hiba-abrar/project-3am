#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <string>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <cmath>
#include <thread>
#include <chrono>
#include <curl/curl.h>
#include "json.hpp"
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using json = nlohmann::json;

// ROAD NODES - Actual intersections and path points on campus roads
map<int, pair<double, double>> roadNodes = {
    // Main circular road nodes
    {1, {24.9330, 67.1100}},   // Near Admin Block
    {2, {24.9328, 67.1110}},   // Near Library
    {3, {24.9325, 67.1120}},   // Mid circular road
    {4, {24.9320, 67.1130}},   // Near fountain
    {5, {24.9315, 67.1140}},   // Mid section
    {6, {24.9310, 67.1145}},   // Near Street 1
    {7, {24.9305, 67.1140}},   // Survey Lab area
    {8, {24.9305, 67.1130}},   // Civil area
    {9, {24.9310, 67.1120}},   // Back to civil
    {10, {24.9315, 67.1110}},  // Near mechanical
    
    // Inner road nodes
    {11, {24.9320, 67.1115}},  // Inner junction 1
    {12, {24.9318, 67.1125}},  // Inner junction 2
    {13, {24.9315, 67.1130}},  // Inner junction 3
    {14, {24.9312, 67.1135}},  // Inner junction 4
    
    // CSIT area nodes
    {15, {24.9314, 67.1140}},  // CSIT entrance junction
    {16, {24.9311, 67.1138}},  // CSIT offices junction
    {17, {24.9309, 67.1142}},  // NPO Circuit junction
    
    // Ground area nodes
    {18, {24.9325, 67.1150}},  // Ground road 1
    {19, {24.9330, 67.1155}},  // Ground road 2
    {20, {24.9324, 67.1165}},  // Basketball area
    {21, {24.9326, 67.1167}},  // Tennis court junction
    
    // Additional path nodes
    {22, {24.9323, 67.1142}},  // DMS Cafeteria junction
    {23, {24.9323, 67.1114}},  // Fire lab junction
    {24, {24.9312, 67.1125}},  // Urban Engineering junction
    {25, {24.9308, 67.1138}},  // Convocation area junction
    {26, {24.930785434981235, 67.11380211636057}},
    {27, {24.9313, 67.1141}}   // CSIT Department - Google Maps coordinates
};

// ROAD CONNECTIONS - Define which nodes are connected by walkable roads
map<int, vector<pair<int, double>>> roadGraph;

// CSIT DEPARTMENT COORDINATES - From Google Maps
const pair<double, double> CSIT_DEPARTMENT_COORDINATES = {24.9313, 67.1141};

// NAMED LOCATIONS mapped to nearest road nodes
// GLOBAL PREDEFINED LOCATIONS
map<pair<double, double>, string> predefined_locations = {
    {{24.9328525, 67.1099341}, "NED University Admin Block"},
    {{24.933102716272494, 67.11100098699812}, "NED University Library"},
    {{24.9312322150391, 67.11249466009184}, "Urban Infrastructure Engineering Department"},
    {{24.93139344975187, 67.11276243453311}, "Civil Engineering Class Rooms"},
    {{24.932268374432034, 67.11347525532474}, "Fire Lab"},
    {{24.931374111201748, 67.11398043939803}, "CSIT Labs"},
    {{24.930842519246948, 67.11395439753406}, "CSIT Department Entrance"},
    {{24.9311094, 67.1135478}, "CSIT Offices"},
    {{24.9349057, 67.1107172}, "Polymer and Petrochemical Engineering"},
    {{24.9338711062338, 67.11102216552638}, "Mosque"},
    {{24.9330082, 67.1099861}, "NED Circular Road"},
    {{24.931090370212964, 67.11450622340925}, "Street 1"},
    {{24.9299167, 67.1156389}, "NED University Main Gate"},
    {{24.9325, 67.1139}, "Ring Street"},
    {{24.9325, 67.1138}, "Road 1"},
    {{24.9322, 67.1142}, "Road 2"},
    {{24.9335, 67.1155}, "NED Ground Road"},
    {{24.933304625311184, 67.11512899716485}, "Ground Road"},
    {{24.932403048156, 67.1165460313555}, "Basketball Court"},
    {{24.93255462400844, 67.11662939578912}, "Tennis Court"},
    {{24.93275629837479, 67.11568797739794}, "Futsal Court"},
    {{24.932141195231143, 67.11631129009899}, "Athletics Track"},
    {{24.932149578523966, 67.116210098944}, "Football Ground"},
    {{24.930802833974123, 67.11522084971281}, "Hockey Ground"},
    {{24.932376114986397, 67.11460777770577}, "Cricket Ground"},
    {{24.9302112, 67.1138518}, "Convocation Ground"},
    {{24.9307957, 67.1131117}, "Civil AV Hall"},
    {{24.93093833913481, 67.11301775226582}, "Civil Lecture Hall"},
    {{24.93190295462124, 67.1126189323467}, "Main Auditorium"},
    {{24.9317433, 67.1128339}, "Fountain Area"},
    {{24.931827918349338, 67.11302423350993}, "Girls Common Room"},
    {{24.932394422988278, 67.1141342937413}, "DMS Cafeteria"},
    {{24.932293372446324, 67.11422188338926}, "Meezan Bank ATM"},
    {{24.933137028815153, 67.1148424432239}, "Girls Gym"},
    {{24.933336990339118, 67.11548757561029}, "Boys Gym"},
    {{24.92896, 67.11352}, "NED Visitor Gate"},
    {{24.92964263413819, 67.11302834820181}, "National Incubation Centre"},
    {{24.93039196072071, 67.11233596592935}, "NED Service Department"},
    {{24.930496161679198, 67.11241634531031}, "SFC Stationary"},
    {{24.930534410454527, 67.11220801681318}, "SFC Canteen"},
    {{24.930696368061145, 67.11234746621153}, "Urban Lawn"},
    {{24.931419456714742, 67.11188542304869}, "Mechanical Lawn"},
    {{24.9314, 67.1121}, "NED Staff Centre"},
    {{24.9315, 67.1115}, "Mech Corner Cafe"},
    {{24.932181120075303, 67.11088250500814}, "NED Medical Centre"},
    {{24.9326, 67.1112}, "STEM Centre"},
    {{24.932413256559116, 67.11083249162313}, "Dean Office"},
    {{24.932895524967876, 67.11051947487461}, "NED White House"},
    {{24.93358, 67.10963}, "Transport Section"},
    {{24.93114915331957, 67.11356200753656}, "Mathematics Department"},
    {{24.93161946024292, 67.1119738417881}, "Mechanical Engineering Department"},
    {{24.934111146559626, 67.11265619152542}, "Environmental Engineering"},
    {{24.931873806908865, 67.11239793449437}, "Electrical Engineering Department"}};

struct Node {
    string name;
    double lat, lon;
};

struct GPSPosition {
    double lat, lon;
    double accuracy;
};

struct NavigationStep {
    string instruction;
    double distance_to_next;
    string direction;
    pair<double, double> location;
    bool is_turn_point;
    string maneuver_type;
};

const double TURN_DETECTION_DISTANCE = 15.0;
const double GPS_UPDATE_INTERVAL = 3.0;

GPSPosition current_device_gps = {0, 0, 0};
bool gps_active = false;

void speak(const string &text);
double haversine(double lat1, double lon1, double lat2, double lon2);
double calculate_bearing(double lat1, double lon1, double lat2, double lon2);
string get_turn_direction(double bearing1, double bearing2);
void get_device_gps_location();
void initialize_road_graph();
vector<NavigationStep> verify_with_google_maps(double start_lat, double start_lon, double end_lat, double end_lon);
size_t WriteCallback(void* contents, size_t size, size_t nmemb, string* userp);

void speak(const string &text) {
    cout << "🎯 " << text << endl;
    this_thread::sleep_for(chrono::milliseconds(1200));
}

double haversine(double lat1, double lon1, double lat2, double lon2) {
    const double R = 6371000;
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    lat1 *= M_PI / 180.0;
    lat2 *= M_PI / 180.0;
    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1) * cos(lat2) * sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;
}

double calculate_bearing(double lat1, double lon1, double lat2, double lon2) {
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    lat1 *= M_PI / 180.0;
    lat2 *= M_PI / 180.0;
    double y = sin(dLon) * cos(lat2);
    double x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon);
    double bearing = atan2(y, x) * 180.0 / M_PI;
    return fmod((bearing + 360.0), 360.0);
}

string get_turn_direction(double bearing1, double bearing2) {
    double diff = bearing2 - bearing1;
    if (diff > 180) diff -= 360;
    if (diff < -180) diff += 360;
    
    if (abs(diff) < 30) return "straight";
    else if (diff > 0) return "right";
    else return "left";
}

void get_device_gps_location() {
    cout << "\n📍 Attempting to get GPS location..." << endl;
    
    cout << "Enter current latitude (or press Enter to use last known): ";
    string lat_input;
    getline(cin, lat_input);
    
    if (!lat_input.empty()) {
        current_device_gps.lat = stod(lat_input);
        
        cout << "Enter current longitude: ";
        string lon_input;
        getline(cin, lon_input);
        current_device_gps.lon = stod(lon_input);
        
        current_device_gps.accuracy = 5.0;
        gps_active = true;
        
        cout << "✅ GPS Location acquired: " << fixed << setprecision(6) 
             << current_device_gps.lat << ", " << current_device_gps.lon << endl;
    }
}

void initialize_road_graph() {
    // Build the road network by connecting nodes that have roads between them
    // Format: roadGraph[fromNode] = {{toNode, distance}, ...}
    
    // Main circular road connections
    roadGraph[1] = {{2, haversine(roadNodes[1].first, roadNodes[1].second, roadNodes[2].first, roadNodes[2].second)},
                    {10, haversine(roadNodes[1].first, roadNodes[1].second, roadNodes[10].first, roadNodes[10].second)}};
    
    roadGraph[2] = {{1, haversine(roadNodes[2].first, roadNodes[2].second, roadNodes[1].first, roadNodes[1].second)},
                    {3, haversine(roadNodes[2].first, roadNodes[2].second, roadNodes[3].first, roadNodes[3].second)},
                    {11, haversine(roadNodes[2].first, roadNodes[2].second, roadNodes[11].first, roadNodes[11].second)}};
    
    roadGraph[3] = {{2, haversine(roadNodes[3].first, roadNodes[3].second, roadNodes[2].first, roadNodes[2].second)},
                    {4, haversine(roadNodes[3].first, roadNodes[3].second, roadNodes[4].first, roadNodes[4].second)},
                    {23, haversine(roadNodes[3].first, roadNodes[3].second, roadNodes[23].first, roadNodes[23].second)}};
    
    roadGraph[4] = {{3, haversine(roadNodes[4].first, roadNodes[4].second, roadNodes[3].first, roadNodes[3].second)},
                    {5, haversine(roadNodes[4].first, roadNodes[4].second, roadNodes[5].first, roadNodes[5].second)},
                    {12, haversine(roadNodes[4].first, roadNodes[4].second, roadNodes[12].first, roadNodes[12].second)}};
    
    roadGraph[5] = {{4, haversine(roadNodes[5].first, roadNodes[5].second, roadNodes[4].first, roadNodes[4].second)},
                    {6, haversine(roadNodes[5].first, roadNodes[5].second, roadNodes[6].first, roadNodes[6].second)},
                    {15, haversine(roadNodes[5].first, roadNodes[5].second, roadNodes[15].first, roadNodes[15].second)},
                    {22, haversine(roadNodes[5].first, roadNodes[5].second, roadNodes[22].first, roadNodes[22].second)}};
    
    roadGraph[6] = {{5, haversine(roadNodes[6].first, roadNodes[6].second, roadNodes[5].first, roadNodes[5].second)},
                    {7, haversine(roadNodes[6].first, roadNodes[6].second, roadNodes[7].first, roadNodes[7].second)},
                    {17, haversine(roadNodes[6].first, roadNodes[6].second, roadNodes[17].first, roadNodes[17].second)}};
    
    roadGraph[7] = {{6, haversine(roadNodes[7].first, roadNodes[7].second, roadNodes[6].first, roadNodes[6].second)},
                    {8, haversine(roadNodes[7].first, roadNodes[7].second, roadNodes[8].first, roadNodes[8].second)}};
    
    roadGraph[8] = {{7, haversine(roadNodes[8].first, roadNodes[8].second, roadNodes[7].first, roadNodes[7].second)},
                    {9, haversine(roadNodes[8].first, roadNodes[8].second, roadNodes[9].first, roadNodes[9].second)},
                    {25, haversine(roadNodes[8].first, roadNodes[8].second, roadNodes[25].first, roadNodes[25].second)}};
    
    roadGraph[9] = {{8, haversine(roadNodes[9].first, roadNodes[9].second, roadNodes[8].first, roadNodes[8].second)},
                    {10, haversine(roadNodes[9].first, roadNodes[9].second, roadNodes[10].first, roadNodes[10].second)},
                    {24, haversine(roadNodes[9].first, roadNodes[9].second, roadNodes[24].first, roadNodes[24].second)}};
    
    roadGraph[10] = {{9, haversine(roadNodes[10].first, roadNodes[10].second, roadNodes[9].first, roadNodes[9].second)},
                     {1, haversine(roadNodes[10].first, roadNodes[10].second, roadNodes[1].first, roadNodes[1].second)}};
    
    // Inner road connections
    roadGraph[11] = {{2, haversine(roadNodes[11].first, roadNodes[11].second, roadNodes[2].first, roadNodes[2].second)},
                     {12, haversine(roadNodes[11].first, roadNodes[11].second, roadNodes[12].first, roadNodes[12].second)}};
    
    roadGraph[12] = {{11, haversine(roadNodes[12].first, roadNodes[12].second, roadNodes[11].first, roadNodes[11].second)},
                     {4, haversine(roadNodes[12].first, roadNodes[12].second, roadNodes[4].first, roadNodes[4].second)},
                     {13, haversine(roadNodes[12].first, roadNodes[12].second, roadNodes[13].first, roadNodes[13].second)}};
    
    roadGraph[13] = {{12, haversine(roadNodes[13].first, roadNodes[13].second, roadNodes[12].first, roadNodes[12].second)},
                     {14, haversine(roadNodes[13].first, roadNodes[13].second, roadNodes[14].first, roadNodes[14].second)}};
    
    roadGraph[14] = {{13, haversine(roadNodes[14].first, roadNodes[14].second, roadNodes[13].first, roadNodes[13].second)},
                     {16, haversine(roadNodes[14].first, roadNodes[14].second, roadNodes[16].first, roadNodes[16].second)}};
    
    // CSIT area connections
    roadGraph[15] = {{5, haversine(roadNodes[15].first, roadNodes[15].second, roadNodes[5].first, roadNodes[5].second)},
                     {16, haversine(roadNodes[15].first, roadNodes[15].second, roadNodes[16].first, roadNodes[16].second)},
                     {22, haversine(roadNodes[15].first, roadNodes[15].second, roadNodes[22].first, roadNodes[22].second)}};
    
    roadGraph[16] = {{15, haversine(roadNodes[16].first, roadNodes[16].second, roadNodes[15].first, roadNodes[15].second)},
                     {14, haversine(roadNodes[16].first, roadNodes[16].second, roadNodes[14].first, roadNodes[14].second)},
                     {17, haversine(roadNodes[16].first, roadNodes[16].second, roadNodes[17].first, roadNodes[17].second)}};
    
    roadGraph[17] = {{16, haversine(roadNodes[17].first, roadNodes[17].second, roadNodes[16].first, roadNodes[16].second)},
                     {6, haversine(roadNodes[17].first, roadNodes[17].second, roadNodes[6].first, roadNodes[6].second)},
                     {25, haversine(roadNodes[17].first, roadNodes[17].second, roadNodes[25].first, roadNodes[25].second)}};
    
    // Ground area connections
    roadGraph[18] = {{23, haversine(roadNodes[18].first, roadNodes[18].second, roadNodes[23].first, roadNodes[23].second)},
                     {19, haversine(roadNodes[18].first, roadNodes[18].second, roadNodes[19].first, roadNodes[19].second)}};
    
    roadGraph[19] = {{18, haversine(roadNodes[19].first, roadNodes[19].second, roadNodes[18].first, roadNodes[18].second)},
                     {20, haversine(roadNodes[19].first, roadNodes[19].second, roadNodes[20].first, roadNodes[20].second)}};
    
    roadGraph[20] = {{19, haversine(roadNodes[20].first, roadNodes[20].second, roadNodes[19].first, roadNodes[19].second)},
                     {21, haversine(roadNodes[20].first, roadNodes[20].second, roadNodes[21].first, roadNodes[21].second)}};
    
    roadGraph[21] = {{20, haversine(roadNodes[21].first, roadNodes[21].second, roadNodes[20].first, roadNodes[20].second)}};
    
    // Additional connections
    roadGraph[22] = {{15, haversine(roadNodes[22].first, roadNodes[22].second, roadNodes[15].first, roadNodes[15].second)},
                     {5, haversine(roadNodes[22].first, roadNodes[22].second, roadNodes[5].first, roadNodes[5].second)}};
    
    roadGraph[23] = {{3, haversine(roadNodes[23].first, roadNodes[23].second, roadNodes[3].first, roadNodes[3].second)},
                     {18, haversine(roadNodes[23].first, roadNodes[23].second, roadNodes[18].first, roadNodes[18].second)},
                     {24, haversine(roadNodes[23].first, roadNodes[23].second, roadNodes[24].first, roadNodes[24].second)}};
    
    roadGraph[24] = {{23, haversine(roadNodes[24].first, roadNodes[24].second, roadNodes[23].first, roadNodes[23].second)},
                     {9, haversine(roadNodes[24].first, roadNodes[24].second, roadNodes[9].first, roadNodes[9].second)}};
    
    roadGraph[25] = {{8, haversine(roadNodes[25].first, roadNodes[25].second, roadNodes[8].first, roadNodes[8].second)},
                     {17, haversine(roadNodes[25].first, roadNodes[25].second, roadNodes[17].first, roadNodes[17].second)}};
}

vector<int> dijkstra_shortest_path(int start_node, int end_node) {
    map<int, double> distances;
    map<int, int> previous;
    set<int> unvisited;
    
    for (const auto& node : roadGraph) {
        distances[node.first] = numeric_limits<double>::infinity();
        unvisited.insert(node.first);
    }
    distances[start_node] = 0;
    
    while (!unvisited.empty()) {
        int current = -1;
        double min_dist = numeric_limits<double>::infinity();
        
        for (int node : unvisited) {
            if (distances[node] < min_dist) {
                min_dist = distances[node];
                current = node;
            }
        }
        
        if (current == -1 || current == end_node) break;
        
        unvisited.erase(current);
        
        for (const auto& neighbor : roadGraph[current]) {
            int next_node = neighbor.first;
            double edge_weight = neighbor.second;
            double alt = distances[current] + edge_weight;
            
            if (alt < distances[next_node]) {
                distances[next_node] = alt;
                previous[next_node] = current;
            }
        }
    }
    
    vector<int> path;
    int current = end_node;
    
    while (current != start_node && previous.find(current) != previous.end()) {
        path.push_back(current);
        current = previous[current];
    }
    path.push_back(start_node);
    reverse(path.begin(), path.end());
    
    return path;
}

vector<NavigationStep> calculate_road_route(const string& start_name, const string& end_name) {
    vector<NavigationStep> route;
    
    if (locationToNode.find(start_name) == locationToNode.end() || 
        locationToNode.find(end_name) == locationToNode.end()) {
        speak("Location not found in road network");
        return route;
    }
    
    int start_node = locationToNode[start_name];
    int end_node = locationToNode[end_name];
    
    vector<int> path = dijkstra_shortest_path(start_node, end_node);
    
    if (path.empty() || path.size() < 2) {
        speak("No route found between locations");
        return route;
    }
    
    speak("Found route through " + to_string(path.size()) + " road nodes");
    
    for (size_t i = 0; i < path.size(); i++) {
        NavigationStep step;
        int node_id = path[i];
        step.location = roadNodes[node_id];
        
        if (i == 0) {
            step.maneuver_type = "depart";
            step.direction = "straight";
            step.is_turn_point = false;
            step.instruction = "Head straight";
        }
        else if (i == path.size() - 1) {
            step.maneuver_type = "arrive";
            step.direction = "straight";
            step.is_turn_point = false;
            step.instruction = "Destination reached";
        }
        else {
            double bearing1 = calculate_bearing(
                roadNodes[path[i-1]].first, roadNodes[path[i-1]].second,
                roadNodes[path[i]].first, roadNodes[path[i]].second
            );
            double bearing2 = calculate_bearing(
                roadNodes[path[i]].first, roadNodes[path[i]].second,
                roadNodes[path[i+1]].first, roadNodes[path[i+1]].second
            );
            
            step.direction = get_turn_direction(bearing1, bearing2);
            step.is_turn_point = (step.direction != "straight");
            
            double distance_to_next = haversine(
                roadNodes[path[i]].first, roadNodes[path[i]].second,
                roadNodes[path[i+1]].first, roadNodes[path[i+1]].second
            );
            
            if (step.is_turn_point) {
                step.instruction = "Turn " + step.direction;
                step.maneuver_type = "turn";
            } else {
                step.instruction = "Walk straight";
                step.maneuver_type = "continue";
            }
        }
        
        if (i < path.size() - 1) {
            step.distance_to_next = haversine(
                roadNodes[path[i]].first, roadNodes[path[i]].second,
                roadNodes[path[i+1]].first, roadNodes[path[i+1]].second
            );
        } else {
            step.distance_to_next = 0;
        }
        
        route.push_back(step);
    }
    
    return route;
}

void visualize_road_route(const vector<NavigationStep>& route, const string& start_name, const string& end_name) {
    cout << "\n🗺  ROAD NODE-BASED ROUTE (Like Google Maps Walking)" << endl;
    cout << "From: " << start_name << endl;
    cout << "To: " << end_name << endl;
    cout << "Path: ";
    
    for (size_t i = 0; i < route.size(); i++) {
        if (i == 0) {
            cout << "🚩START";
        }
        else if (i == route.size() - 1) {
            cout << "🎯END";
        }
        else {
            if (route[i].is_turn_point) {
                if (route[i].direction == "left") {
                    cout << "⬅LEFT";
                } else if (route[i].direction == "right") {
                    cout << "➡RIGHT";
                }
            } else {
                cout << "⬆STRAIGHT";
            }
        }
        
        if (i < route.size() - 1) {
            cout << " → ";
        }
    }
    cout << endl;
    
    cout << "\n🔄 DETAILED TURN INSTRUCTIONS:" << endl;
    for (size_t i = 0; i < route.size(); i++) {
        if (i == 0) {
            cout << "🚩 " << route[i].instruction << endl;
        }
        else if (i == route.size() - 1) {
            cout << "🎯 " << route[i].instruction << endl;
        }
        else {
            cout << "➡ " << route[i].instruction;
            if (route[i].distance_to_next > 0) {
                cout << " (" << int(route[i].distance_to_next) << "m)";
            }
            cout << endl;
        }
    }
    cout << endl;
}

void real_time_road_navigation(const vector<NavigationStep>& route, const string& dest_name) {
    speak("Starting navigation with REAL GPS tracking on roads");
    speak("Your device location will be used for navigation");
    
    int current_segment = 0;
    bool destination_reached = false;
    vector<bool> turn_announced(route.size(), false);
    
    cout << "\n🚶 STARTING REAL-TIME GPS NAVIGATION" << endl;
    cout << "📍 Using device GPS location..." << endl;
    
    speak(route[0].instruction);
    
    while (!destination_reached && current_segment < route.size() - 1) {
        get_device_gps_location();
        
        if (!gps_active || current_device_gps.lat == 0) {
            cout << "⚠  Waiting for GPS signal..." << endl;
            this_thread::sleep_for(chrono::seconds(int(GPS_UPDATE_INTERVAL)));
            continue;
        }
        
        double target_lat = route[current_segment + 1].location.first;
        double target_lon = route[current_segment + 1].location.second;
        
        double distance_to_next = haversine(current_device_gps.lat, current_device_gps.lon, 
                                           target_lat, target_lon);
        
        for (size_t i = current_segment + 1; i < route.size() - 1; i++) {
            if (route[i].is_turn_point && !turn_announced[i]) {
                double distance_to_turn = 0;
                for (int j = current_segment; j < i; j++) {
                    distance_to_turn += route[j].distance_to_next;
                }
                distance_to_turn -= distance_to_next;
                
                if (distance_to_turn <= TURN_DETECTION_DISTANCE) {
                    speak(route[i].instruction);
                    turn_announced[i] = true;
                    cout << "🔄 TURN ANNOUNCED: " << route[i].instruction << endl;
                }
            }
        }
        
        if (current_segment < route.size() - 1 && !route[current_segment].is_turn_point) {
            if (distance_to_next > 30 && fmod(distance_to_next, 30.0) < 2.0) {
                string progress = "Continue straight for " + to_string(int(distance_to_next)) + " meters";
                speak(progress);
            }
        }
        
        if (distance_to_next < 10.0) {
            if (current_segment == route.size() - 2) {
                speak("You have arrived at " + dest_name);
                speak("Navigation complete!");
                destination_reached = true;
                cout << "🎯 DESTINATION REACHED: " << dest_name << endl;
                break;
            } else {
                current_segment++;
                cout << "📍 Reached navigation point " << current_segment << endl;
                
                if (current_segment < route.size() - 1 && !route[current_segment].is_turn_point) {
                    speak("Continue straight");
                }
            }
        }
        
        cout << "📡 DEVICE GPS: " << fixed << setprecision(6) 
             << current_device_gps.lat << ", " << current_device_gps.lon 
             << " | To next: " << int(distance_to_next) << "m"
             << " | Segment: " << current_segment + 1 << "/" << route.size() 
             << " | Accuracy: ±" << int(current_device_gps.accuracy) << "m" << endl;
        
        this_thread::sleep_for(chrono::milliseconds(int(GPS_UPDATE_INTERVAL * 1000)));
    }
    
    if (!destination_reached) {
        speak("Navigation ended");
    }
}

void wait_for_space() {
    cout << "Press SPACE then ENTER to continue...";
    string input;
    getline(cin, input);
}

int get_numeric_input(int min_val, int max_val) {
    int choice;
    while (true) {
        if (cin >> choice) {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            if (choice >= min_val && choice <= max_val) {
                return choice;
            }
        } else {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        speak("Invalid input. Please enter a number between " + to_string(min_val) + " and " + to_string(max_val));
    }
}

string get_location_from_user(const string &prompt) {
    speak(prompt);
    cout << "Available locations:" << endl;
    int count = 1;
    for (const auto& loc : locationToNode) {
        cout << count++ << ". " << loc.first << endl;
    }
    
    cout << "\nType location name: ";
    string location_name;
    getline(cin, location_name);
    
    for (const auto& loc : locationToNode) {
        if (loc.first.find(location_name) != string::npos || 
            location_name.find(loc.first) != string::npos) {
            return loc.first;
        }
    }
    
    return "";
}

void start_road_navigation() {
    speak("NED University Navigation System - Road Node-Based Routing");
    speak("Like Google Maps walking mode - using actual road intersections");
    
    string start_name = get_location_from_user("Enter starting location");
    if (start_name.empty()) {
        speak("Invalid starting location");
        return;
    }
    
    string dest_name = get_location_from_user("Enter destination location");
    if (dest_name.empty()) {
        speak("Invalid destination location");
        return;
    }
    
    speak("Calculating optimal route from " + start_name + " to " + dest_name);
    
    vector<NavigationStep> route = calculate_road_route(start_name, dest_name);
    
    if (route.empty()) {
        speak("Could not calculate route. Please try different locations.");
        return;
    }
    
    visualize_road_route(route, start_name, dest_name);
    
    double total_distance = 0;
    for (const auto& step : route) {
        total_distance += step.distance_to_next;
    }
    
    speak("Route calculated. Total distance: " + to_string(int(total_distance)) + " meters");
    speak("Starting real-time GPS navigation in 3 seconds...");
    speak("Make sure GPS is enabled on your device");
    this_thread::sleep_for(chrono::seconds(3));
    
    real_time_road_navigation(route, dest_name);
}

void list_locations_menu() {
    cout << "\n📍 AVAILABLE LOCATIONS:" << endl;
    int count = 1;
    for (const auto& loc : locationToNode) {
        cout << "  " << count++ << ". " << loc.first << " (Node " << loc.second << ")" << endl;
    }
}

int main() {
    initialize_road_graph();
    
    speak("NED University Navigation System Started");
    speak("Voice-guided navigation with REAL GPS tracking using road nodes");
    speak("Like Google Maps walking mode");
    
    while (true) {
        cout << "\n=== MAIN MENU ===" << endl;
        cout << "1. List Locations" << endl;
        cout << "2. Start Navigation" << endl;
        cout << "3. Exit" << endl;
        cout << "Choice: ";
        
        int choice = get_numeric_input(1, 3);
        
        switch (choice) {
            case 1:
                list_locations_menu();
                wait_for_space();
                break;
            case 2:
                start_road_navigation();
                wait_for_space();
                break;
            case 3:
                speak("Thank you for using NED Navigation. Goodbye!");
                return 0;
        }
    }
    
    return 0;
}