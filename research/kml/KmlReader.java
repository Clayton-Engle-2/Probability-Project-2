package research.kml;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import com.vividsolutions.jts.geom.Coordinate;

import research.model.Municipality;

public class KmLReader {

    private String kmlFile;

    public KmLReader(String kmlFile) {
        this.kmlFile = kmlFile;
    }

    public Map<String, Municipality> readKml() {
        Map<String, Municipality> municipalities = new HashMap<String, Municipality>();
        try {
            XMLInputFactory inputFactory = XMLInputFactory.newInstance();
            FileInputStream inputStream = new FileInputStream(kmlFile);
            XMLStreamReader reader = inputFactory.createXMLStreamReader(inputStream);

            Municipality municipality = null;

            while (reader.hasNext()) {
                int eventType = reader.next();

                if (eventType == XMLStreamReader.START_ELEMENT) {
                    String elementName = reader.getLocalName();

                    if (elementName.equalsIgnoreCase("Municipal_Boundaries_of_NJ")) {
                        municipality = new Municipality();
                    } else if (municipality != null) {
                        if (elementName.equalsIgnoreCase("OBJECTID")) {
                            municipality.setObjectID(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("MUN")) {
                            municipality.setName(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("COUNTY")) {
                            municipality.setCounty(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("MUN_LABEL")) {
                            municipality.setLabel(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("MUN_TYPE")) {
                            municipality.setType(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("GNIS_NAME")) {
                            municipality.setGnisName(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("GNIS")) {
                            municipality.setGnis(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("SSN")) {
                            municipality.setSsn(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("MUN_CODE")) {
                            municipality.setMunCode(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("CENSUS2010")) {
                            municipality.setCensus2010(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("ACRES")) {
                            municipality.setAcres(Double.parseDouble(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("SQ_MILES")) {
                            municipality.setSqMiles(Double.parseDouble(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("POP2010")) {
                            municipality.setPop2010(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("POP2000")) {
                            municipality.setPop2000(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("POP1990")) {
                            municipality.setPop1990(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("POP1980")) {
                            municipality.setPop1980(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("POPDEN2010")) {
                            municipality.setPopden2010(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("POPDEN2000")) {
                            municipality.setPopden2000(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("POPDEN1990")) {
                            municipality.setPopden1990(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("POPDEN1980")) {
                            municipality.setPopden1980(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("Shape_Length")) {
                            municipality.setShapeLength(Double.parseDouble(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("Shape_Area")) {
                            municipality.setShapeArea(Double.parseDouble(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("CENSUS2020")) {
                            municipality.setCensus2020(reader.getElementText());
                        } else if (elementName.equalsIgnoreCase("POP2020")) {
                            municipality.setPop2020(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("POPDEN2020")) {
                            municipality.setPopden2020(Integer.parseInt(reader.getElementText()));
                        } else if (elementName.equalsIgnoreCase("coordinates")) {
                            String coordinatesText = reader.getElementText();
                            String[] coordinatePairs = coordinatesText.trim().split("\\s+");
                            List<Coordinate> coordinates = new ArrayList<>();

                            for (String coordinatePair : coordinatePairs) {
                                String[] coordinateValues = coordinatePair.split(",");
                                double longitude = Double.parseDouble(coordinateValues[0]);
                                double latitude = Double.parseDouble(coordinateValues[1]);
                                Coordinate coordinate = new Coordinate(latitude, longitude);
                                coordinates.add(coordinate);
                            }
                            municipality.setBoundaryCoordinates(coordinates);
                        }
                    }

                } else if (eventType == XMLStreamReader.END_ELEMENT) {
                	if (reader.getLocalName().equalsIgnoreCase("Municipal_Boundaries_of_NJ")) {
                	    municipalities.put(municipality.getCounty() + municipality.getName(), municipality);
                	    municipality = null;
                	}
                }
            }

            reader.close();
        } catch (FileNotFoundException | XMLStreamException e) {
            e.printStackTrace();
        }

        return municipalities;
    }
}

