package research.csv;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import research.model.District;
import research.model.Municipality;

public class CSVReaderUtility {

    public static void readCSVAndUpdateMunicipalities(String csvFile, HashMap<String, Municipality> municipalities, District[] districts) {
        try (FileReader fileReader = new FileReader(csvFile); BufferedReader bufferedReader = new BufferedReader(fileReader)) {
            String line;
            int currentDistrict = -1;
            boolean header = true;

            while ((line = bufferedReader.readLine()) != null) {
                if (header) {
                    header = false;
                    continue;
                }

                String[] values = line.split(",");

                if (values[0].startsWith("District")) {
                    currentDistrict = Integer.parseInt(values[0].split(" ")[1]) - 1;
                    continue;
                }

                String county = values[0].trim();
                String name = values[1].trim();
                double demVotes = Double.parseDouble(values[2]);
                double repVotes = Double.parseDouble(values[3]);
                double othVotes = Double.parseDouble(values[4]);

                for (int i = 5; i < values.length; i++) {
                    othVotes += Double.parseDouble(values[i]);
                }

                String municipalityKey = county + name;
                Municipality municipality = municipalities.get(municipalityKey);

                if (municipality != null) {
                    municipality.setDemVotes(demVotes);
                    municipality.setRepVotes(repVotes);
                    municipality.setOthVotes(othVotes);
                    districts[currentDistrict].addTown(municipality);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
