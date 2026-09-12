PRAGMA foreign_keys = ON;

DROP VIEW  IF EXISTS v_analysis;
DROP VIEW  IF EXISTS v_sample_annotated;
DROP VIEW  IF EXISTS v_sample_frequencies;
DROP TABLE IF EXISTS sample_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS treatment_courses;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS cell_populations;
DROP TABLE IF EXISTS treatments;
DROP TABLE IF EXISTS conditions;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
  project_id   INTEGER PRIMARY KEY,
  project_code TEXT NOT NULL UNIQUE
);

CREATE TABLE conditions (
  condition_id   INTEGER PRIMARY KEY,
  condition_name TEXT NOT NULL UNIQUE
);

CREATE TABLE treatments (
  treatment_id   INTEGER PRIMARY KEY,
  treatment_name TEXT NOT NULL UNIQUE
);

CREATE TABLE subjects (
  subject_id   INTEGER PRIMARY KEY,
  subject_code TEXT    NOT NULL UNIQUE,
  project_id   INTEGER NOT NULL REFERENCES projects(project_id),
  condition_id INTEGER NOT NULL REFERENCES conditions(condition_id),
  age          INTEGER CHECK (age IS NULL OR age BETWEEN 0 AND 120),
  sex          TEXT    NOT NULL CHECK (sex IN ('M','F'))
);

-- One row per subject per drug. Response is a property of the course,
-- not of an individual aliquot. NULL response = untreated or not yet assessed.
CREATE TABLE treatment_courses (
  course_id    INTEGER PRIMARY KEY,
  subject_id   INTEGER NOT NULL REFERENCES subjects(subject_id),
  treatment_id INTEGER NOT NULL REFERENCES treatments(treatment_id),
  response     TEXT CHECK (response IS NULL OR response IN ('yes','no')),
  UNIQUE (subject_id, treatment_id)
);

CREATE TABLE samples (
  sample_id                 INTEGER PRIMARY KEY,
  sample_code               TEXT    NOT NULL UNIQUE,
  course_id                 INTEGER NOT NULL REFERENCES treatment_courses(course_id),
  sample_type               TEXT    NOT NULL,
  time_from_treatment_start INTEGER NOT NULL,
  UNIQUE (course_id, sample_type, time_from_treatment_start)
);

CREATE TABLE cell_populations (
  population_id   INTEGER PRIMARY KEY,
  population_name TEXT NOT NULL UNIQUE
);

-- Long format. One row per sample per population.
CREATE TABLE sample_counts (
  sample_id     INTEGER NOT NULL REFERENCES samples(sample_id) ON DELETE CASCADE,
  population_id INTEGER NOT NULL REFERENCES cell_populations(population_id),
  count         INTEGER NOT NULL CHECK (count >= 0),
  PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX idx_subjects_project    ON subjects(project_id);
CREATE INDEX idx_subjects_condition  ON subjects(condition_id);
CREATE INDEX idx_subjects_sex        ON subjects(sex);
CREATE INDEX idx_courses_subject     ON treatment_courses(subject_id);
CREATE INDEX idx_courses_treatment   ON treatment_courses(treatment_id);
CREATE INDEX idx_courses_response    ON treatment_courses(response);
CREATE INDEX idx_samples_course      ON samples(course_id);
CREATE INDEX idx_samples_type_time   ON samples(sample_type, time_from_treatment_start);
CREATE INDEX idx_counts_population   ON sample_counts(population_id, sample_id);

-- Part 2 deliverable. Column names are exactly as the assignment specifies.
CREATE VIEW v_sample_frequencies AS
SELECT
  s.sample_code AS sample,
  SUM(sc.count) OVER (PARTITION BY sc.sample_id)                            AS total_count,
  cp.population_name                                                        AS population,
  sc.count                                                                  AS count,
  ROUND(100.0 * sc.count / SUM(sc.count) OVER (PARTITION BY sc.sample_id), 4) AS percentage
FROM sample_counts sc
JOIN samples          s  ON s.sample_id      = sc.sample_id
JOIN cell_populations cp ON cp.population_id = sc.population_id;

CREATE VIEW v_sample_annotated AS
SELECT
  s.sample_id,
  s.sample_code  AS sample,
  s.sample_type,
  s.time_from_treatment_start,
  sub.subject_code AS subject,
  sub.age,
  sub.sex,
  p.project_code   AS project,
  c.condition_name AS condition,
  t.treatment_name AS treatment,
  tc.response
FROM samples s
JOIN treatment_courses tc  ON tc.course_id    = s.course_id
JOIN subjects          sub ON sub.subject_id  = tc.subject_id
JOIN projects          p   ON p.project_id    = sub.project_id
JOIN conditions        c   ON c.condition_id  = sub.condition_id
JOIN treatments        t   ON t.treatment_id  = tc.treatment_id;

-- Single entry point for Parts 3 and 4 and the dashboard.
CREATE VIEW v_analysis AS
SELECT a.*, f.population, f.count, f.total_count, f.percentage
FROM v_sample_annotated a
JOIN v_sample_frequencies f ON f.sample = a.sample;