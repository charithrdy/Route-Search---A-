#!/usr/bin/env python
import string
from collections import defaultdict
from Queue import PriorityQueue
from math import radians, cos, sin, asin, sqrt
import sys


def solve_ids(start_city, end_city, algo, cost):
    for i in range(10000):
        fringe = []
        distance_so_far = 0
        time_so_far = 0
        route_so_far = start_city
        fringe.append((start_city, route_so_far, 0, distance_so_far, time_so_far))
        visited = defaultdict(list)
        while len(fringe):
            if fringe[len(fringe) - 1][2] <= i:
                (state, route_so_far, depth_so_far, distance_so_far, time_so_far) = fringe.pop()
                if state == end_city:
                    if cost == "segments":
                        return "yes " + str(distance_so_far) + " " + str(time_so_far) + " " + route_so_far
                    else:
                        return "no " + str(distance_so_far) + " " + str(time_so_far) + " " + route_so_far
                if not visited[state]:
                    visited[state] = True
                    if depth_so_far + 1 <= i:
                        for city in succOfCity[state]:
                            fringe.append((city[0], route_so_far + " " + city[0], depth_so_far + 1,
                                           distance_so_far + city[1], time_so_far + float(city[1]) / float(city[2])))
    return False


def solve_dfs(start_city, end_city, algo, cost):
    fringe = []
    route_so_far = start_city
    distance_so_far = 0
    time_so_far = 0
    fringe.append((start_city, route_so_far, distance_so_far, time_so_far))
    visited = defaultdict(list)
    while len(fringe):
        (state, route_so_far, distance_so_far, time_so_far) = fringe.pop()
        if state == end_city:
            return "no " + str(distance_so_far) + " " + str(time_so_far) + " " + route_so_far
        if not visited[state]:
            visited[state] = True
            for city in succOfCity[state]:
                fringe.append((city[0], route_so_far + " " + city[0], distance_so_far + city[1],
                               time_so_far + float(city[1]) / float(city[2])))
    return False


def solve_bfs(start_city, end_city, algo, cost):
    fringe = []
    route_so_far = start_city
    distance_so_far = 0
    time_so_far = 0
    fringe.append((start_city, route_so_far, distance_so_far, time_so_far))
    visited = defaultdict(list)
    while len(fringe):
        (state, route_so_far, distance_so_far, time_so_far) = fringe.pop(0)
        if state == end_city:
            if cost == "segments":
                return "yes " + str(distance_so_far) + " " + str(time_so_far) + " " + route_so_far
            else:
                return "no " + str(distance_so_far) + " " + str(time_so_far) + " " + route_so_far
        if not visited[state]:
            visited[state] = True
            for city in succOfCity[state]:
                fringe.append((city[0], route_so_far + " " + city[0], distance_so_far + city[1],
                               time_so_far + float(city[1]) / float(city[2])))
    return False


def solve_ucs(start_city, end_city, algo, cost):
    fringe = PriorityQueue()
    route_so_far = start_city
    cost_so_far = 0
    distance_so_far = 0
    time_so_far = 0
    fringe.put((0, (start_city, route_so_far, cost_so_far, distance_so_far, time_so_far)))
    visited = defaultdict(list)
    visited[start_city] = True
    while fringe.qsize() > 0:
        (state, route_so_far, cost_so_far, distance_so_far, time_so_far) = fringe.get()[1]
        visited[state] = True
        if state == end_city:
            return "yes " + str(distance_so_far) + " " + str(time_so_far) + " " + route_so_far
        else:
            for city in succOfCity[state]:
                if not visited[city[0]]:
                    if cost == "distance":
                        current_cost = city[1]
                        fringe.put((current_cost + cost_so_far, (
                        city[0], route_so_far + " " + city[0], current_cost + cost_so_far, distance_so_far + city[1],
                        time_so_far + float(city[1]) / float(city[2]))))
                    elif cost == "time":
                        current_cost = float(city[1]) / float(city[2])
                        fringe.put((current_cost + cost_so_far, (
                        city[0], route_so_far + " " + city[0], current_cost + cost_so_far, distance_so_far + city[1],
                        time_so_far + float(city[1]) / float(city[2]))))
                    elif cost == "segments":
                        current_cost = 1
                        fringe.put((current_cost + cost_so_far, (
                        city[0], route_so_far + " " + city[0], current_cost + cost_so_far, distance_so_far + city[1],
                        time_so_far + float(city[1]) / float(city[2]))))
    return False


def calcDisplacement(fromCity, toCity):
    """
        Haversine method to Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)
    """
    # Referred from the below link.
    # https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
    global prev_state
    if latlon[fromCity]:
        lon1 = latlon[fromCity][0]
        lat1 = latlon[fromCity][1]
    else:
        # If lat, lon fof cities missing

        if fromCity == "" or len(citiesInState[fromCity.split(",")[1]]) == 0:
            # If nearest city with lon, lat found in it's succesors and unable to calc state average because of invalid name or has no cities in the state with lat, lon data
            lon1 = LonsofState[prev_state] / len(citiesInState[prev_state])
            lat1 = LatsofState[prev_state] / len(citiesInState[prev_state])
        else:
            # If nearest city with lon, lat found in it's succesors, calclate state avg
            lon1 = LonsofState[fromCity.split(",")[1]] / len(citiesInState[fromCity.split(",")[1]])
            lat1 = LatsofState[fromCity.split(",")[1]] / len(citiesInState[fromCity.split(",")[1]])
    if latlon[toCity]:
        lon2 = latlon[toCity][0]
        lat2 = latlon[toCity][1]
    else:
        lon2 = LonsofState[toCity.split(",")[1]] / len(citiesInState[toCity.split(",")[1]])
        lat2 = LatsofState[toCity.split(",")[1]] / len(citiesInState[toCity.split(",")[1]])

    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # global prev_state
    if fromCity != "":
        prev_state = fromCity.split(",")[1]
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 3956  # Radius of earth in kilometers. Use 3956 for miles
    return c * r


def findClosestCity(city):
    min = 10000
    nearestCity = ""
    roadDist = 0
    for nearbyCities in succOfCity[city]:
        if latlon[nearbyCities[0]]:
            roadDist = nearbyCities[1]
            if roadDist < min:
                min = roadDist
                nearestCity = nearbyCities[0]
        else:
            continue
    return nearestCity


def heuristic_astar(city):
    if latlon[city] and latlon[end_city]:
        return calcDisplacement(city, end_city)
    elif latlon[city] and not latlon[end_city]:
        return calcDisplacement(city, findClosestCity(end_city))
    elif not latlon[city] and latlon[end_city]:
        return calcDisplacement(findClosestCity(city), end_city)
    else:
        return calcDisplacement(findClosestCity(city), findClosestCity(end_city))


def solve_astar(start_city, end_city, algo, cost):
    fringe = PriorityQueue()
    route_so_far = start_city
    cost_so_far = 0
    distance_so_far = 0
    time_so_far = 0
    fringe.put((0, (start_city, route_so_far, cost_so_far, distance_so_far, time_so_far)))
    visited = defaultdict(list)
    visited[start_city] = True
    while fringe.qsize() > 0:
        (state, route_so_far, cost_so_far, distance_so_far, time_so_far) = fringe.get()[1]
        visited[state] = True
        if state == end_city:
            return "yes " + str(distance_so_far) + " " + str(time_so_far) + " " + route_so_far
        else:
            for city in succOfCity[state]:
                if not visited[city[0]]:
                    if cost == "distance":
                        current_cost = city[1]
                        fringe.put((current_cost + cost_so_far + heuristic_astar(city[0]), (
                        city[0], route_so_far + " " + city[0], current_cost + cost_so_far, distance_so_far + city[1],
                        time_so_far + float(city[1]) / float(city[2]))))
                    elif cost == "time":
                        current_cost = float(city[1]) / float(city[2])
                        fringe.put((current_cost + cost_so_far + float(heuristic_astar(city[0])) / float(avg_speed), (
                        city[0], route_so_far + " " + city[0], current_cost + cost_so_far, distance_so_far + city[1],
                        time_so_far + float(city[1]) / float(city[2]))))
                    elif cost == "segments":
                        current_cost = 1
                        fringe.put((current_cost + cost_so_far + heuristic_astar(city[0]) / float(avg_dist), (
                        city[0], route_so_far + " " + city[0], current_cost + cost_so_far, distance_so_far + city[1],
                        time_so_far + float(city[1]) / float(city[2]))))
    return False


succOfCity = defaultdict(list)
speedsofHighway = defaultdict(list)
sumOfSpeeds = 0
numOfRoutes = 0
sumOfDist = 0

with open('road-segments.txt', 'r') as file:
    for line in file:
        if not len(line.split()) == 4 and int(line.split()[3]) != 0:
            succOfCity[line.split()[0]].append(
                (line.split()[1], int(line.split()[2]), int(line.split()[3]), line.split()[4]))
            succOfCity[line.split()[1]].append(
                (line.split()[0], int(line.split()[2]), int(line.split()[3]), line.split()[4]))
            speedsofHighway[line.split()[4]].append(int(line.split()[3]))
        numOfRoutes = numOfRoutes + 1

with open('road-segments.txt', 'r') as file:
    for line in file:
        sumOfDist = sumOfDist + int(line.split()[2])
        if len(line.split()) == 4 or int(line.split()[3]) == 0:
            if not speedsofHighway[line.split()[3]]:
                speedsofHighway[line.split()[3]].append(int(30))
                sumOfSpeeds = sumOfSpeeds + 30
            else:
                sumOfSpeeds = sumOfSpeeds + int(
                    sum(speedsofHighway[line.split()[3]]) / len(speedsofHighway[line.split()[3]]))
            succOfCity[line.split()[0]].append((line.split()[1], int(line.split()[2]), int(
                sum(speedsofHighway[line.split()[3]]) / len(speedsofHighway[line.split()[3]])), line.split()[3]))
            succOfCity[line.split()[1]].append((line.split()[0], int(line.split()[2]), int(
                sum(speedsofHighway[line.split()[3]]) / len(speedsofHighway[line.split()[3]])), line.split()[3]))
        else:
            sumOfSpeeds = sumOfSpeeds + int(line.split()[3])

avg_speed = float(sumOfSpeeds) / float(numOfRoutes)
avg_dist = float(sumOfDist) / float(numOfRoutes)
latlon = defaultdict(list)
citiesInState = defaultdict(list)
LonsofState = defaultdict(float)
LatsofState = defaultdict(float)

with open('city-gps.txt', 'r') as file:
    for line in file:
        citiesInState[(line.split()[0]).split(",")[1]].append(line.split()[0])
        if not LonsofState[(line.split()[0]).split(",")[1]]:
            LonsofState[(line.split()[0]).split(",")[1]] = 0
        if not LatsofState[(line.split()[0]).split(",")[1]]:
            LatsofState[(line.split()[0]).split(",")[1]] = 0
        LonsofState[(line.split()[0]).split(",")[1]] = LonsofState[(line.split()[0]).split(",")[1]] + float(
            line.split()[1])
        LatsofState[(line.split()[0]).split(",")[1]] = LatsofState[(line.split()[0]).split(",")[1]] + float(
            line.split()[2])
        latlon[line.split()[0]].append((float(line.split()[1])))
        latlon[line.split()[0]].append((float(line.split()[2])))

start_city = sys.argv[1]
end_city = sys.argv[2]
algo = sys.argv[3]
cost = sys.argv[4]
prev_state = start_city.split(",")[1]
if algo == "bfs":
    print(solve_bfs(start_city, end_city, algo, cost))
if algo == "uniform":
    print(solve_ucs(start_city, end_city, algo, cost))
if algo == "dfs":
    print(solve_dfs(start_city, end_city, algo, cost))
if algo == "ids":
    print(solve_ids(start_city, end_city, algo, cost))
if algo == "astar":
    print(solve_astar(start_city, end_city, algo, cost))