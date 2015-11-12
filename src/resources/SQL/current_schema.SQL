-- MySQL dump 10.13  Distrib 5.5.46, for debian-linux-gnu (x86_64)
--
-- Host: 127.0.0.1    Database: RNA_refinery
-- ------------------------------------------------------
-- Server version	5.5.46-0ubuntu0.12.04.2

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `annotation_cdhit_genomic`
--

DROP TABLE IF EXISTS `annotation_cdhit_genomic`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_cdhit_genomic` (
  `query` varchar(256) NOT NULL,
  `hit` varchar(256) NOT NULL,
  `affiliation` varchar(256) DEFAULT NULL,
  `desc` varchar(1024) NOT NULL,
  `identity` float DEFAULT NULL,
  `bits` int(11) DEFAULT NULL,
  `qlen` int(11) DEFAULT NULL,
  `alen` int(11) DEFAULT NULL,
  `gaps` int(11) DEFAULT NULL,
  `evalue` float DEFAULT NULL,
  `run_id` int(11) DEFAULT NULL,
  FULLTEXT KEY `cdhit_genomic_query` (`query`),
  FULLTEXT KEY `cdhit_genomic_hit` (`hit`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotation_cdhit_pep`
--

DROP TABLE IF EXISTS `annotation_cdhit_pep`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_cdhit_pep` (
  `query` varchar(256) NOT NULL,
  `hit` varchar(256) NOT NULL,
  `affiliation` varchar(256) DEFAULT NULL,
  `desc` varchar(1024) NOT NULL,
  `identity` float DEFAULT NULL,
  `bits` int(11) DEFAULT NULL,
  `qlen` int(11) DEFAULT NULL,
  `alen` int(11) DEFAULT NULL,
  `gaps` int(11) DEFAULT NULL,
  `evalue` float DEFAULT NULL,
  `run_id` int(11) DEFAULT NULL,
  FULLTEXT KEY `cdhit_pep_query` (`query`),
  FULLTEXT KEY `cdhit_pep_hit` (`hit`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotation_cdhit_rna`
--

DROP TABLE IF EXISTS `annotation_cdhit_rna`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_cdhit_rna` (
  `query` varchar(256) NOT NULL,
  `hit` varchar(256) NOT NULL,
  `affiliation` varchar(256) DEFAULT NULL,
  `desc` varchar(1024) NOT NULL,
  `identity` float DEFAULT NULL,
  `bits` int(11) DEFAULT NULL,
  `qlen` int(11) DEFAULT NULL,
  `alen` int(11) DEFAULT NULL,
  `gaps` int(11) DEFAULT NULL,
  `evalue` float DEFAULT NULL,
  `run_id` int(11) DEFAULT NULL,
  FULLTEXT KEY `cdhit_rna_query` (`query`),
  FULLTEXT KEY `cdhit_rna_hit` (`hit`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotation_cdna`
--

DROP TABLE IF EXISTS `annotation_cdna`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_cdna` (
  `qName` varchar(256) NOT NULL,
  `score` int(11) DEFAULT NULL,
  `qStart` int(11) DEFAULT NULL,
  `qEnd` int(11) DEFAULT NULL,
  `qStrand` char(1) DEFAULT NULL,
  `qLength` int(11) DEFAULT NULL,
  `qCoverage` float DEFAULT NULL,
  `identity` float DEFAULT NULL,
  `similarity` float DEFAULT NULL,
  `tName` varchar(256) NOT NULL,
  `tStart` int(11) DEFAULT NULL,
  `tEnd` int(11) DEFAULT NULL,
  `tStrand` char(1) DEFAULT NULL,
  `tLength` int(11) DEFAULT NULL,
  `tCoverage` float DEFAULT NULL,
  `qNumInsert` int(11) DEFAULT NULL,
  `qBaseInsert` int(11) DEFAULT NULL,
  `tNumInsert` int(11) DEFAULT NULL,
  `tBaseInsert` int(11) DEFAULT NULL,
  `run_id` int(11) DEFAULT NULL,
  FULLTEXT KEY `cdna_qName` (`qName`),
  FULLTEXT KEY `cdna_tName` (`tName`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotation_cds`
--

DROP TABLE IF EXISTS `annotation_cds`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_cds` (
  `qName` varchar(256) NOT NULL,
  `score` int(11) DEFAULT NULL,
  `qStart` int(11) DEFAULT NULL,
  `qEnd` int(11) DEFAULT NULL,
  `qStrand` char(1) DEFAULT NULL,
  `qLength` int(11) DEFAULT NULL,
  `qCoverage` float DEFAULT NULL,
  `identity` float DEFAULT NULL,
  `similarity` float DEFAULT NULL,
  `tName` varchar(256) NOT NULL,
  `tStart` int(11) DEFAULT NULL,
  `tEnd` int(11) DEFAULT NULL,
  `tStrand` char(1) DEFAULT NULL,
  `tLength` int(11) DEFAULT NULL,
  `tCoverage` float DEFAULT NULL,
  `qNumInsert` int(11) DEFAULT NULL,
  `qBaseInsert` int(11) DEFAULT NULL,
  `tNumInsert` int(11) DEFAULT NULL,
  `tBaseInsert` int(11) DEFAULT NULL,
  `run_id` int(11) DEFAULT NULL,
  FULLTEXT KEY `cds_qName` (`qName`),
  FULLTEXT KEY `cds_tName` (`tName`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotation_genome`
--

DROP TABLE IF EXISTS `annotation_genome`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_genome` (
  `qName` varchar(256) NOT NULL,
  `score` int(11) DEFAULT NULL,
  `qStart` int(11) DEFAULT NULL,
  `qEnd` int(11) DEFAULT NULL,
  `qStrand` char(1) DEFAULT NULL,
  `qLength` int(11) DEFAULT NULL,
  `qCoverage` float DEFAULT NULL,
  `identity` float DEFAULT NULL,
  `similarity` float DEFAULT NULL,
  `tName` varchar(256) NOT NULL,
  `tStart` int(11) DEFAULT NULL,
  `tEnd` int(11) DEFAULT NULL,
  `tStrand` char(1) DEFAULT NULL,
  `tLength` int(11) DEFAULT NULL,
  `tCoverage` float DEFAULT NULL,
  `qNumInsert` int(11) DEFAULT NULL,
  `qBaseInsert` int(11) DEFAULT NULL,
  `tNumInsert` int(11) DEFAULT NULL,
  `tBaseInsert` int(11) DEFAULT NULL,
  `run_id` int(11) DEFAULT NULL,
  FULLTEXT KEY `genome_qName` (`qName`),
  FULLTEXT KEY `genome_tName` (`tName`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotation_refseq_dna`
--

DROP TABLE IF EXISTS `annotation_refseq_dna`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_refseq_dna` (
  `query` varchar(256) NOT NULL,
  `hit` varchar(256) NOT NULL,
  `desc` varchar(1024) NOT NULL,
  `identity` float DEFAULT NULL,
  `bits` int(11) DEFAULT NULL,
  `qlen` int(11) DEFAULT NULL,
  `alen` int(11) DEFAULT NULL,
  `gaps` int(11) DEFAULT NULL,
  `evalue` float DEFAULT NULL,
  `run_id` int(11) DEFAULT NULL,
  FULLTEXT KEY `ref_dna_query` (`query`),
  FULLTEXT KEY `ref_dna_hit` (`hit`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotation_refseq_pep`
--

DROP TABLE IF EXISTS `annotation_refseq_pep`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_refseq_pep` (
  `query` varchar(256) NOT NULL,
  `hit` varchar(256) NOT NULL,
  `desc` varchar(1024) NOT NULL,
  `identity` float DEFAULT NULL,
  `bits` int(11) DEFAULT NULL,
  `qlen` int(11) DEFAULT NULL,
  `alen` int(11) DEFAULT NULL,
  `gaps` int(11) DEFAULT NULL,
  `evalue` float DEFAULT NULL,
  `run_id` int(11) DEFAULT NULL,
  FULLTEXT KEY `ref_pep_query` (`query`),
  FULLTEXT KEY `ref_pep_hit` (`hit`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `centroids`
--

DROP TABLE IF EXISTS `centroids`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `centroids` (
  `centroid` varchar(256) NOT NULL,
  `centroid_run_id` int(11) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `clusters`
--

DROP TABLE IF EXISTS `clusters`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `clusters` (
  `centroid` varchar(256) NOT NULL,
  `centroid_run_id` int(11) DEFAULT NULL,
  `member` varchar(256) NOT NULL,
  `member_run_id` int(11) DEFAULT NULL,
  `similarity` float DEFAULT NULL,
  KEY `centroid_idx` (`centroid_run_id`),
  KEY `member_idx` (`member_run_id`),
  FULLTEXT KEY `centr_name` (`centroid`),
  FULLTEXT KEY `seq_name` (`member`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `raw_contigs`
--

DROP TABLE IF EXISTS `raw_contigs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `raw_contigs` (
  `run_id` int(11) DEFAULT NULL,
  `name` varchar(256) NOT NULL,
  `sequence` longtext NOT NULL,
  `is_excluded` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`name`),
  UNIQUE KEY `seq_name` (`name`),
  KEY `runid_key` (`run_id`),
  FULLTEXT KEY `ctg_name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sra_runs`
--

DROP TABLE IF EXISTS `sra_runs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sra_runs` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `run_accession` varchar(32) DEFAULT NULL,
  `raw_reads` int(11) DEFAULT NULL,
  `contigs` int(11) DEFAULT NULL,
  `tissue` varchar(256) DEFAULT NULL,
  `condition` varchar(256) DEFAULT NULL,
  `dev_stage` varchar(256) DEFAULT NULL,
  `meta_data` varchar(256) DEFAULT NULL,
  `processed` datetime DEFAULT '0000-00-00 00:00:00',
  PRIMARY KEY (`id`),
  UNIQUE KEY `run_name` (`run_accession`)
) ENGINE=MyISAM AUTO_INCREMENT=312 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2015-11-12 11:09:30
