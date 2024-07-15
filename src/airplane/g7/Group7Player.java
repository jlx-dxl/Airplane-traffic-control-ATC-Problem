package airplane.g7;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Random;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.IntStream;

import org.apache.log4j.Logger;
import airplane.sim.Plane;
import airplane.sim.Player;

public class Group7Player extends Player {

    private final Logger logger = Logger.getLogger(this.getClass()); // for logging
    private Random random = new Random();
    private ArrayList<Airline> airlines = new ArrayList<>();
    private ArrayList<Airline> bestAirlines = new ArrayList<>();

    private final double alpha = 10.0; // time weight
    private final double beta = 2.0; // power weight
    private final double gamma = 10.0; // delay weight

    private final double coolingRate = 1-1e-3;
    private final boolean ifBend = true;

    private final double curvature = 0.3;
    private final int maxPossibleDelay = 1500;

    private final double distanceThreshold = 6.0;
    private final double collisionPenalty = 1e8;

    private double best_cost = 0.0;
    private double bestTimeCost = 0.0;
    private double bestPowerCost = 0.0;
    private double bestDelayCost = 0.0;

    private double temperature = 1000;

    private class Airline {
        private int departureTime;
        private int arrivalTime;
        private ArrayList<Point2D.Double> path;
        private ArrayList<Double> bearings;

        public Airline() {
            this.departureTime = 0;
            this.arrivalTime = 0;
            this.path = new ArrayList<Point2D.Double>();
            this.bearings = new ArrayList<Double>();
        }

        public Airline(int departureTime, ArrayList<Point2D.Double> path, ArrayList<Double> bearings) {
            this.departureTime = departureTime;
            this.arrivalTime = departureTime + path.size();
            this.path = path;
            this.bearings = bearings;
        }

        public boolean isFlying(int round) {
            return (round >= departureTime && round < arrivalTime);
        }

        @Override
        public String toString() {
            return String.format("\nDeparture Time: %d, Arrival Time: %d, Path Size: %s, Bearings Size: %s;", departureTime, arrivalTime, path.size(), bearings.size());
        }
    }

    @Override
    public String getName() {
        return "Group7Player";
    }

    @Override
    public void startNewGame(ArrayList<Plane> planes) {
        logger.info("Start Solving ... ");

        // 初始化Airlines，每个Airline对应一个Plane，记录一条直接从起点到终点的路径和航向以及起飞和到达时间
        for (Plane plane : planes) {
            ArrayList<Point2D.Double> path = new ArrayList<>();
            ArrayList<Double> bearings = new ArrayList<>();

            Point2D.Double start = new Point2D.Double(plane.getLocation().getX(), plane.getLocation().getY());
            Point2D.Double end = plane.getDestination();
            double bearing = calculateBearing(start, end);

            Point2D.Double current = start;
            while (!current.equals(end)) {
                path.add(current);
                bearings.add(bearing);
                current = moveTowards(current, bearing);
                if (current.distance(end) < 0.5) {
                    current = end;
                } else if (current.distance(end) > 100) {
                    logger.error("Distance Error!");
                    break;
                }
            }
            path.add(end);
            bearings.add(bearing);

            airlines.add(new Airline(plane.getDepartureTime(), path, bearings));
        }


        // 开始模拟退火

        double cost = 0.0;

        ArrayList<Airline> newAirlines = new ArrayList<>();

        best_cost = calculateCost(planes, airlines);

        while (temperature > 1) {

            // 随机延后未起飞飞机的起飞时间
            newAirlines = randomlyDelayAirlines(airlines);

            // 随机弯曲曲线，更新path，bearings和arrivalTime

            if (ifBend) {
                newAirlines = randomlyBendAirlines(newAirlines, planes);
            }

            // 计算cost
            cost = calculateCost(planes, newAirlines);

            // 根据温度决定是否接收解
            if (acceptanceProbability(best_cost, cost, temperature) > random.nextDouble()) {
                bestAirlines = newAirlines;
                best_cost = cost;
            }

            // 温度下降
            temperature *= coolingRate;

        }
        logger.info(String.format("Best Time Cost: %f, Best Power Cost: %f, Best Delay Cost: %f, Best Total Cost: %f", bestTimeCost, bestPowerCost, bestDelayCost, best_cost));
        logger.info("Best Airlines: " + bestAirlines.toString());
        logger.info("Finish Solving ...");
    }

    @Override
    public double[] updatePlanes(ArrayList<Plane> planes, int round, double[] bearings) {
        for (int i = 0; i < planes.size(); i++) {
            Plane plane = planes.get(i);
            Airline airline = bestAirlines.get(i);

            if (round < airline.departureTime) {
                bearings[i] = -1; // Not yet departed
            } else if (airline.departureTime<=round && bearings[i] != -2) {
                bearings[i] = airline.bearings.get(round - airline.departureTime); // Update to the optimal bearing
            } else {
                bearings[i] = -2; // Arrived at destination
            }
        }
        return bearings;
    }

    private Point2D.Double moveTowards(Point2D.Double point, double bearing) {
        double radialBearing = Math.toRadians(bearing - 90);
        double newX = point.getX() + Math.cos(radialBearing);
        double newY = point.getY() + Math.sin(radialBearing);
        return new Point2D.Double(newX, newY);
    }

    private double calculateCost(ArrayList<Plane>planes, ArrayList<Airline>airlines) {
        double timeCost = 0.0;
        double powerCost = 0.0;
        double delayCost = 0.0;
        double totalCost = 0.0;

        for (int i = 0; i < planes.size(); i++) {
            Plane plane = planes.get(i);
            Airline airline = airlines.get(i);

            // 计算时间代价
            if (airline.arrivalTime - 1 > timeCost) {
                timeCost = airline.arrivalTime - 1;
            }

            // 计算能量代价
            powerCost += airline.path.size() - 1;

            // 计算延迟代价
            if (airline.departureTime > plane.getDepartureTime()) {
                delayCost += (airline.departureTime - plane.getDepartureTime());
            }
        }

        totalCost = alpha * timeCost + beta * powerCost + gamma * delayCost;

        double distance = 0.0;

        // i代表回合数，
        for (int i = 0; i < (int) timeCost; i++) {
            for (int j = 0; j < planes.size(); j++) {
                if (!airlines.get(j).isFlying(i)) {
                    continue;
                } else {
                    for (int k = 0; k < planes.size(); k++) {
                        if (!airlines.get(k).isFlying(i)) {
                            continue;
                        } else {
                            if (j == k) {
                                    continue;
                            } else {
                                distance = airlines.get(j).path.get(i-airlines.get(j).departureTime).distance(airlines.get(k).path.get(i-airlines.get(k).departureTime));
                                if (distance < distanceThreshold) {
                                    totalCost += collisionPenalty;
                                }
                            }
                        }
                    }
                }
            }
        }

        logger.info(String.format("Temperature: %f, Time Cost: %f, Power Cost: %f, Delay Cost: %f, Total Cost: %f", temperature, timeCost, powerCost, delayCost, totalCost));

        if (totalCost < best_cost) {
            bestTimeCost = timeCost;
            bestPowerCost = powerCost;
            bestDelayCost = delayCost;
        }

        return totalCost;
    }

    // 高温阶段:
    // 在算法的初期，温度较高。此时，即使新的解比当前解差，算法也有较高的概率接受这个较差的解。
    // 这使得算法能够进行较大范围的搜索，探索更多的解空间，从而避免陷入局部最优。
    // 低温阶段:
    // 随着算法的进行，温度逐渐降低。此时，只有当新的解显著优于当前解时，算法才会接受这个新的解。
    // 这使得算法在搜索的后期更倾向于细致地优化当前解，逐步收敛到全局最优解。
    private double acceptanceProbability(double currentCost, double newCost, double temperature) {
        if (newCost < currentCost) { // 如果新代价更低，直接接受
            return 1.0;
        }
        return Math.exp((currentCost - newCost) / temperature); // 否则根据温度计算接受概率
    }

    private ArrayList<Airline> randomlyDelayAirlines(ArrayList<Airline> airlines) {
        ArrayList<Airline> delayedAirlines = new ArrayList<>();
        for (Airline airline : airlines) {
            Airline delayedAirline = new Airline(airline.departureTime + random.nextInt(maxPossibleDelay), airline.path, airline.bearings);
            delayedAirlines.add(delayedAirline);
        }
        return delayedAirlines;
    }

    private ArrayList<Airline> randomlyBendAirlines(ArrayList<Airline> airlines, ArrayList<Plane> planes) {
        for (int i = 0; i < planes.size(); i++) {
            Plane plane = planes.get(i);
            Airline airline = airlines.get(i);
            airline.path = generateArcPoints(plane.getLocation(), plane.getDestination(), curvature*(random.nextDouble()*2.0-1.0));
            airline.arrivalTime = airline.departureTime + airline.path.size();
            airline.bearings = generateBearings(airline.path);
        }
        return airlines;
    }

    private ArrayList<Point2D.Double> generateArcPoints(Point2D.Double start, Point2D.Double end, double curvature) {
        ArrayList<Point2D.Double> points = new ArrayList<>();

        double x1 = start.getX();
        double y1 = start.getY();
        double x2 = end.getX();
        double y2 = end.getY();

        // Calculate the midpoint
        double midX = (x1 + x2) / 2;
        double midY = (y1 + y2) / 2;

        // Calculate the control point for the quadratic Bezier curve
        double controlX = midX + (y1 - y2) * curvature;
        double controlY = midY + (x2 - x1) * curvature;

        // Generate points along the curve using the Bezier formula
        double tStep = 1.0 / calculateBezierCurveLength(start, new Point2D.Double(controlX, controlY), end, 1000); // step size based on distance between points
        for (double t = 0; t < 1.0; t += tStep) {
            double x = (1 - t) * (1 - t) * x1 + 2 * (1 - t) * t * controlX + t * t * x2;
            double y = (1 - t) * (1 - t) * y1 + 2 * (1 - t) * t * controlY + t * t * y2;
            points.add(new Point2D.Double(x, y));
        }

        points.add(end);

        return points;
    }

    private ArrayList<Double> generateBearings(ArrayList<Point2D.Double> points) {
        ArrayList<Double> bearings = new ArrayList<>();

        // Calculate bearings between consecutive points
        for (int i = 0; i < points.size(); i++) {
            if (i == points.size() - 1) {
                bearings.add(bearings.get(i - 1));
            } else {
                bearings.add(calculateBearing(points.get(i), points.get(i + 1)));
            }
        }
        return bearings;
    }

    public static double calculateBezierCurveLength(Point2D.Double start, Point2D.Double control, Point2D.Double end, int numSegments) {
        DoubleUnaryOperator dxdt = t -> 2 * (1 - t) * (control.x - start.x) + 2 * t * (end.x - control.x);
        DoubleUnaryOperator dydt = t -> 2 * (1 - t) * (control.y - start.y) + 2 * t * (end.y - control.y);

        return IntStream.range(0, numSegments)
                .mapToDouble(i -> {
                    double t1 = (double) i / numSegments;
                    double t2 = (double) (i + 1) / numSegments;
                    double dx1 = dxdt.applyAsDouble(t1);
                    double dy1 = dydt.applyAsDouble(t1);
                    double dx2 = dxdt.applyAsDouble(t2);
                    double dy2 = dydt.applyAsDouble(t2);
                    return 0.5 * (Math.sqrt(dx1 * dx1 + dy1 * dy1) + Math.sqrt(dx2 * dx2 + dy2 * dy2)) * (t2 - t1);
                })
                .sum();
    }
}
