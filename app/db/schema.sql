
CREATE TABLE user (
  id TEXT PRIMARY KEY,
  name TEXT,
  email TEXT UNIQUE,
  is_guest INTEGER DEFAULT 0,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

/* table separator */

CREATE TABLE task (
  id TEXT NOT NULL,
  user_id TEXT NOT NULL,
  task_type TEXT NOT NULL,
  status TEXT NOT NULL,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  /*updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,*/
  PRIMARY KEY (id)
  FOREIGN KEY (user_id) REFERENCES user(id)
);

/* todo: fix commented out line (for updated_at) to work with sqlite accepted syntax
    (see: https://stackoverflow.com/questions/6578439/on-update-current-timestamp-with-sqlite)
*/
